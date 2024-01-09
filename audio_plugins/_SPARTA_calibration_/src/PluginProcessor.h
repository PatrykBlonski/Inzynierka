/*
 ==============================================================================
 
 This file is part of SPARTA; a suite of spatial audio plug-ins.
 Copyright (c) 2018 - Leo McCormack.
 
 SPARTA is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 SPARTA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with SPARTA.  If not, see <http://www.gnu.org/licenses/>.
 
 ==============================================================================
*/
#ifndef PLUGINPROCESSOR_H_INCLUDED
#define PLUGINPROCESSOR_H_INCLUDED

#include <JuceHeader.h>
#include "calibration.h"
#include "powermap.h"
#define CONFIGURATIONHELPER_ENABLE_GENERICLAYOUT_METHODS 1
#include "../../resources/ConfigurationHelper.h"
#include <thread>
#include "kiss_fft.h"
#define BUILD_VER_SUFFIX ""  /* String to be added before the version name on the GUI (e.g. beta, alpha etc..) */
#ifndef MIN
# define MIN(a,b) (( (a) < (b) ) ? (a) : (b))
#endif
#ifndef MAX
# define MAX(a,b) (( (a) > (b) ) ? (a) : (b))
#endif

typedef enum _TIMERS{
    TIMER_PROCESSING_RELATED = 1,
    TIMER_GUI_RELATED
}TIMERS;

/* Parameter tags: for the default VST GUI */
enum { 
    k_roomCoeff,
    k_numInputs,
    k_numOutputs,
    
    k_NumOfParameters
};




class PluginProcessor  : public AudioProcessor,
                         public MultiTimer,
                         public VSTCallbackHandler
{
public:
    AudioBuffer<float> recordingBuffer;  // A buffer to store the microphone input.
    AudioBuffer<float> ImpulseBuffer;
    AudioBuffer<float> SweepBuffer;
    AudioBuffer<float> TempBuffer;
    //AudioProcessorValueTreeState parameters;
    /* Get functions */
    void* getFXHandle() { return hPan; }
    void* getFXHandlePmap() { return hPm; }
	bool getIsPlaying() { return isPlaying; }
    int getCurrentBlockSize(){ return nHostBlockSize; }
    int getCurrentNumInputs(){ return nNumInputs; }
    int getCurrentNumOutputs(){ return nNumOutputs; }
    
    /* For refreshing window during automation */
    bool refreshWindow;
    void setRefreshWindow(bool newState) { refreshWindow = newState; }
    bool getRefreshWindow() { return refreshWindow; }
    
    /* JSON */
    void saveConfigurationToFile (File destination, int srcOrLs);
    void loadConfiguration (const File& presetFile, int srcOrLs);
    void setLastDir(File newLastDir){ lastDir = newLastDir; }
    File getLastDir() {return lastDir;}
    
    /* VST CanDo */
    pointer_sized_int handleVstManufacturerSpecific (int32 /*index*/, pointer_sized_int /*value*/, void* /*ptr*/, float /*opt*/) override { return 0; }
    pointer_sized_int handleVstPluginCanDo (int32 /*index*/, pointer_sized_int /*value*/, void* ptr, float /*opt*/) override{
        auto text = (const char*) ptr;
        auto matches = [=](const char* s) { return strcmp (text, s) == 0; };
        if (matches ("wantsChannelCountNotifications"))
            return 1;
        return 0;
    }
    
private:
    void* hPan;           /* calibration handle */
    void* hPm;          /* powermap calibration handle */
    int nNumInputs;       /* current number of input channels */
    int nNumOutputs;      /* current number of output channels */
    int nSampleRate;      /* current host sample rate */
    int nHostBlockSize;   /* typical host block size to expect, in samples */
    File lastDir;
    ValueTree elements {"SourcesOrLoudspeakers"}; 
    
    //new variables
    bool calibrating = false;  // To check if calibration is active
    double phase = 0.0;        // Phase for our sine wave generation
    double frequency = 20.0;   // Starting frequency
    const double startFrequency = 20.0;
    const double endFrequency = 24000.0;
    const double duration = 10.0;   // Duration in seconds
    double timeElapsed = 0.0;  // Keep track of how long the sweep has been playing
    bool isRecording = false;
    dsp::Convolution convolution1;
    dsp::Convolution convolution2;
    dsp::Convolution convolution3;
    dsp::Convolution convolution4;
    int latency = 0;


    void timerCallback(int timerID) override
    {
        switch(timerID){
            case TIMER_PROCESSING_RELATED:
                /* reinitialise codec if needed */
                if(calibration_getCodecStatus(hPan) == CODEC_STATUS_NOT_INITIALISED){
                    try{
                        std::thread threadInit(calibration_initCodec, hPan);
                        threadInit.detach();
                    } catch (const std::exception& exception) {
                        std::cout << "Could not create thread" << exception.what() << std::endl;
                    }
                }
                if (powermap_getCodecStatus(hPm) == CODEC_STATUS_NOT_INITIALISED) {
                    try {
                        std::thread threadInit(powermap_initCodec, hPm);
                        threadInit.detach();
                    }
                    catch (const std::exception& exception) {
                        std::cout << "Could not create thread" << exception.what() << std::endl;
                    }
                }
                break;
            case TIMER_GUI_RELATED:
                /* handled in PluginEditor */
                break;
        }
    }
    
    /***************************************************************************\
                                    JUCE Functions
    \***************************************************************************/
public:
    PluginProcessor();
    ~PluginProcessor();
    bool isPlaying;
    int currentRecordingPosition = 0;
    int loudspeakerNumber = 0;
    bool playAll = false;
    double computeSweepFrequency(double time);
    void PluginProcessor::startCalibration();
    void PluginProcessor::endCalibration();
    void PluginProcessor::generateSine(const double deltaT, AudioBuffer<float>& buffer);
    void PluginProcessor::saveBufferToWav(AudioBuffer<float>& buffer, const String& path);
    AudioBuffer<float> loadImpulseResponse(const String& filePath);
    void loadImpulseResponse();
    void PluginProcessor::distanceCalculation(AudioBuffer<float>& sweep, AudioBuffer<float>& input, int loudNum);
    void PluginProcessor::generateSineSweep(float sampleRate, juce::AudioBuffer<float>& sweepBuffer);
    std::vector<float> dists;
    void PluginProcessor::next();
    


    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;
    void processBlock (AudioSampleBuffer& buffer, MidiBuffer& midiMessages) override;
    AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;
    const String getName() const override;
    int getNumParameters() override;
    float getParameter (int index) override;
    void setParameter (int index, float newValue) override;
    const String getParameterName (int index) override;
    const String getParameterText (int index) override;
    const String getInputChannelName (int channelIndex) const override;
    const String getOutputChannelName (int channelIndex) const override;
    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool silenceInProducesSilenceOut() const override;
    double getTailLengthSeconds() const override;
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram(int index) override;
    const String getProgramName(int index) override;
    bool isInputChannelStereoPair (int index) const override;
    bool isOutputChannelStereoPair(int index) const override;
    void changeProgramName(int index, const String& newName) override;
    void getStateInformation (MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;


private:
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (PluginProcessor)
};

#endif  // PLUGINPROCESSOR_H_INCLUDED
