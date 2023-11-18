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
#include "PluginProcessor.h"
#include "PluginEditor.h"

AudioProcessorValueTreeState::ParameterLayout createParameterLayout()
{
    AudioProcessorValueTreeState::ParameterLayout layout;

    for (int i = 1; i < 9; ++i)
        layout.add(std::make_unique<AudioParameterInt>(String(i), String(i), 0, i, 0));

    return layout;
}

PluginProcessor::PluginProcessor() : 
	AudioProcessor(BusesProperties()
		.withInput("Input", AudioChannelSet::discreteChannels(64), true)
	    .withOutput("Output", AudioChannelSet::discreteChannels(64), true))
{
    AudioProcessorValueTreeState parameters(*this, nullptr, "PARAMETERS", createParameterLayout());
    refreshWindow = true;
    panner_create(&hPan);
    startTimer(TIMER_PROCESSING_RELATED, 80); 
}

PluginProcessor::~PluginProcessor()
{
	panner_destroy(&hPan);
}

void PluginProcessor::startCalibration() {
    calibrating = true;
    phase = 0.0;
    timeElapsed = 0.0;
    frequency = startFrequency;
    isRecording = true;
    currentRecordingPosition = 0;
    loudspeakerNumber = 0;
}

void PluginProcessor::endCalibration() {
    loudspeakerNumber++;
    if (loudspeakerNumber >= panner_getNumLoudspeakers(hPan)) {
        calibrating = 0;
        isRecording = false;
    }
    else {
        phase = 0.0;
        timeElapsed = 0.0;
        frequency = startFrequency;
        currentRecordingPosition = 0;
    }
}

double PluginProcessor::computeSweepFrequency(double time) {
    return jmap(time, 0.0, duration * 5, startFrequency, endFrequency);        
}

void PluginProcessor::setParameter (int index, float newValue)
{
    /* standard parameters */
    if(index < k_NumOfParameters){
        switch (index) {
            case k_roomCoeff:  panner_setDTT(hPan, newValue); break;
            case k_numInputs:  panner_setNumSources(hPan, (int)(newValue*(float)(MAX_NUM_INPUTS)+0.5)); break;
            case k_numOutputs: panner_setNumLoudspeakers(hPan, (int)(newValue*(float)(MAX_NUM_OUTPUTS)+0.5)); break;
        }
    }
    /* source direction parameters */
    else if(index<2*MAX_NUM_INPUTS+k_NumOfParameters){
        index-=k_NumOfParameters;
        float newValueScaled;
        if (!(index % 2)){
            newValueScaled = (newValue - 0.5f)*360.0f;
            if (newValueScaled != panner_getSourceAzi_deg(hPan, index/2)){
                panner_setSourceAzi_deg(hPan, index/2, newValueScaled);
                refreshWindow = true;
            }
        }
        else{
            newValueScaled = (newValue - 0.5f)*180.0f;
            if (newValueScaled != panner_getSourceElev_deg(hPan, index/2)){
                panner_setSourceElev_deg(hPan, index/2, newValueScaled);
                refreshWindow = true;
            }
        }
    }
    /* loudspeaker direction parameters */
    else {
        index -= (k_NumOfParameters + 2 * MAX_NUM_INPUTS);
        float newValueScaled;
        int loudspeakerIndex = index / 3;

        if (index % 3 == 0) {
            // Azimuth
            newValueScaled = (newValue - 0.5f) * 360.0f;
            if (newValueScaled != panner_getLoudspeakerAzi_deg(hPan, loudspeakerIndex)) {
                panner_setLoudspeakerAzi_deg(hPan, loudspeakerIndex, newValueScaled);
                refreshWindow = true;
            }
        }
        else if (index % 3 == 1) {
            // Elevation
            newValueScaled = (newValue - 0.5f) * 180.0f;
            if (newValueScaled != panner_getLoudspeakerElev_deg(hPan, loudspeakerIndex)) {
                panner_setLoudspeakerElev_deg(hPan, loudspeakerIndex, newValueScaled);
                refreshWindow = true;
            }
        }
        else {
            // Distance
            // Assuming a maximum distance of MAX_DISTANCE
            newValueScaled = newValue * 10;
            if (newValueScaled != panner_getLoudspeakerDist_deg(hPan, loudspeakerIndex)) {
                panner_setLoudspeakerDist_deg(hPan, loudspeakerIndex, newValueScaled);
                refreshWindow = true;
            }
        }
    }
}

void PluginProcessor::setCurrentProgram (int /*index*/)
{
}

float PluginProcessor::getParameter (int index)
{
    /* standard parameters */
    if(index < k_NumOfParameters){
        switch (index) {
            case k_roomCoeff:  return (float)panner_getDTT(hPan);
            case k_numInputs:  return (float)(panner_getNumSources(hPan))/(float)(MAX_NUM_INPUTS);
            case k_numOutputs: return (float)(panner_getNumLoudspeakers(hPan))/(float)(MAX_NUM_OUTPUTS);
            default: return 0.0f;
        }
    }
    /* source direction parameters */
    else if(index<2*MAX_NUM_INPUTS+k_NumOfParameters){
        index-=k_NumOfParameters;
        if (!(index % 2))
            return (panner_getSourceAzi_deg(hPan, index/2)/360.0f) + 0.5f;
        else
            return (panner_getSourceElev_deg(hPan, (index-1)/2)/180.0f) + 0.5f;
    }
    /* loudspeaker direction parameters */
    else{
        index -= (k_NumOfParameters+2*MAX_NUM_INPUTS);
        if (!(index % 2))
            return (panner_getLoudspeakerAzi_deg(hPan, index/2)/360.0f) + 0.5f;
        else
            return (panner_getLoudspeakerElev_deg(hPan, (index-1)/2)/180.0f) + 0.5f;
    }
}

int PluginProcessor::getNumParameters()
{
	return k_NumOfParameters + 2*MAX_NUM_INPUTS + 2*MAX_NUM_OUTPUTS;
}

const String PluginProcessor::getName() const
{
    return JucePlugin_Name;
}

const String PluginProcessor::getParameterName (int index)
{
    /* standard parameters */
    if(index < k_NumOfParameters){
        switch (index) {
            case k_roomCoeff:  return "roomCoeff";
            case k_numInputs:  return "num_sources";
            case k_numOutputs: return "num_loudspeakers";
            default: return "NULL";
        }
    }
    /* source direction parameters */
    else if(index<2*MAX_NUM_INPUTS+k_NumOfParameters){
        index-=k_NumOfParameters;
        if (!(index % 2))
            return TRANS("SrcAzim_") + String(index/2);
        else
            return TRANS("SrcElev_") + String((index-1)/2);
    }
    /* loudspeaker direction parameters */
    else{
        index -= (k_NumOfParameters+2*MAX_NUM_INPUTS);
        if (!(index % 2))
            return TRANS("LsAzim_") + String(index/2);
        else
            return TRANS("LsElev_") + String((index-1)/2);
    }
}

const String PluginProcessor::getParameterText(int index)
{
    /* standard parameters */
    if(index < k_NumOfParameters){
        switch (index) {
            case k_roomCoeff:  return String(panner_getDTT(hPan));
            case k_numInputs:  return String(panner_getNumSources(hPan));
            case k_numOutputs: return String(panner_getNumLoudspeakers(hPan));
            default: return "NULL";
        }
    }
    /* source direction parameters */
    else if(index<2*MAX_NUM_INPUTS+k_NumOfParameters){
        index -= k_NumOfParameters;
        if (!(index % 2))
            return String(panner_getSourceAzi_deg(hPan, index/2));
        else
            return String(panner_getSourceElev_deg(hPan, (index-1)/2));
    }
    /* loudspeaker direction parameters */
    else{
        index -= (k_NumOfParameters+2*MAX_NUM_INPUTS);
        if (!(index % 2))
            return String(panner_getLoudspeakerAzi_deg(hPan, index/2));
        else
            return String(panner_getLoudspeakerElev_deg(hPan, (index-1)/2));
    }
}

const String PluginProcessor::getInputChannelName (int channelIndex) const
{
    return String (channelIndex + 1);
}

const String PluginProcessor::getOutputChannelName (int channelIndex) const
{
    return String (channelIndex + 1);
}

double PluginProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int PluginProcessor::getNumPrograms()
{
    return 0;
}

int PluginProcessor::getCurrentProgram()
{
    return 0;
}

const String PluginProcessor::getProgramName (int /*index*/)
{
    return String();
}


bool PluginProcessor::isInputChannelStereoPair (int /*index*/) const
{
    return true;
}

bool PluginProcessor::isOutputChannelStereoPair (int /*index*/) const
{
    return true;
}


bool PluginProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool PluginProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool PluginProcessor::silenceInProducesSilenceOut() const
{
    return false;
}

void PluginProcessor::changeProgramName (int /*index*/, const String& /*newName*/)
{
}

void PluginProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    nHostBlockSize = samplesPerBlock;
    nNumInputs =  getTotalNumInputChannels();
    nNumOutputs = getTotalNumOutputChannels();
    nSampleRate = (int)(sampleRate + 0.5);
    isPlaying = false;

	panner_init(hPan, nSampleRate);
    AudioProcessor::setLatencySamples(panner_getProcessingDelay());
}

void PluginProcessor::releaseResources()
{
    isPlaying = false;
}

void PluginProcessor::generateSine(const double deltaT, AudioBuffer<float>& buffer) {
    int numberOfInputs = panner_getNumSources(hPan);
    const int numSamples = buffer.getNumSamples();
    for (int sample = 0; sample < numSamples; ++sample) {
        // Compute the instantaneous frequency for the sweep
        frequency = computeSweepFrequency(timeElapsed);

        // Compute the sample of the sine sweep
        float sineSample = std::sin(phase);

        // Write this sample to all output channels (you can change this if needed)
        buffer.setSample(loudspeakerNumber + numberOfInputs, sample, sineSample);

        // Increment phase and wrap if it exceeds 2*PI
        phase += (2.0 * double_Pi * frequency) * deltaT;
        if (phase > 2.0 * double_Pi) {
            phase -= 2.0 * double_Pi;
        }

        // Update time elapsed
        timeElapsed += deltaT;

        // Stop the sweep if the duration has been reached
        if (timeElapsed >= duration) {
            endCalibration();
            break;
        }
    }
}

void PluginProcessor::processBlock(AudioBuffer<float>& buffer, MidiBuffer& midiMessages) {
    // Get some important details
    const int totalNumInputChannels = getTotalNumInputChannels();
    const int totalNumOutputChannels = getTotalNumOutputChannels();
    const int numSamples = buffer.getNumSamples();
    const double sampleRate = getSampleRate();
    const double deltaT = 1.0 / sampleRate;  // Time increment per sample
    int numberOfInputs = panner_getNumSources(hPan);
    recordingBuffer.setSize(numberOfInputs, 48000 * 5);

    // Clear any unused channels
    for (int i = totalNumInputChannels; i < totalNumOutputChannels; ++i) {
        buffer.clear(i, 0, numSamples);
    }



    if (calibrating) {
        generateSine(deltaT, buffer);
    }

    if (isRecording && currentRecordingPosition + buffer.getNumSamples() <= recordingBuffer.getNumSamples()) {
        for (int channel = 0; channel < numberOfInputs; ++channel) {
            recordingBuffer.copyFrom(channel, currentRecordingPosition, buffer, channel, 0, buffer.getNumSamples());
        }
        currentRecordingPosition += buffer.getNumSamples();
    }
}

void PluginProcessor::saveBufferToWav()
{
    // File path
    juce::File outputFile = juce::File::getSpecialLocation(juce::File::userDesktopDirectory).getChildFile("recorded_audio.wav");

    // Create an AudioFormatManager and register the WAV format
    juce::AudioFormatManager formatManager;
    formatManager.registerBasicFormats();  // This registers the WAV format among others

    // Find the WAV format
    juce::AudioFormat* wavFormat = formatManager.findFormatForFileExtension("wav");
    if (wavFormat == nullptr)
    {
        // Error handling: WAV format not found
        return;
    }

    // Create a writer for the specified file and format
    std::unique_ptr<juce::AudioFormatWriter> writer(wavFormat->createWriterFor(
        new juce::FileOutputStream(outputFile),
        getSampleRate(),
        recordingBuffer.getNumChannels(),
        16,  // Bit depth, e.g., 16 bits
        {},  // No metadata
        0    // Use the default quality
    ));

    if (writer.get() != nullptr)
    {
        // Write the entire buffer to the file
        writer->writeFromAudioSampleBuffer(recordingBuffer, 0, recordingBuffer.getNumSamples());
    }
}

//==============================================================================
bool PluginProcessor::hasEditor() const
{
    return true; 
}

AudioProcessorEditor* PluginProcessor::createEditor()
{
    return new PluginEditor (this);
}

//==============================================================================
void PluginProcessor::getStateInformation (MemoryBlock& destData)
{
    XmlElement xml("PANNERPLUGINSETTINGS");
    for(int i=0; i<panner_getMaxNumSources(); i++){
        xml.setAttribute("SourceAziDeg" + String(i), panner_getSourceAzi_deg(hPan,i));
        xml.setAttribute("SourceElevDeg" + String(i), panner_getSourceElev_deg(hPan,i));
    }
    xml.setAttribute("nSources", panner_getNumSources(hPan));
    xml.setAttribute("DTT", panner_getDTT(hPan));
    
    xml.setAttribute("JSONFilePath", lastDir.getFullPathName());
    
    for(int i=0; i<panner_getMaxNumLoudspeakers(); i++){
        xml.setAttribute("LoudspeakerAziDeg" + String(i), panner_getLoudspeakerAzi_deg(hPan,i));
        xml.setAttribute("LoudspeakerElevDeg" + String(i), panner_getLoudspeakerElev_deg(hPan,i));
    }
    xml.setAttribute("nLoudspeakers", panner_getNumLoudspeakers(hPan));

    copyXmlToBinary(xml, destData);
}

void PluginProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    std::unique_ptr<XmlElement> xmlState(getXmlFromBinary(data, sizeInBytes));
    
    if (xmlState != nullptr) {
        if (xmlState->hasTagName("PANNERPLUGINSETTINGS")) {
            for(int i=0; i<panner_getMaxNumSources(); i++){
                if(xmlState->hasAttribute("SourceAziDeg" + String(i)))
                    panner_setSourceAzi_deg(hPan, i, (float)xmlState->getDoubleAttribute("SourceAziDeg" + String(i), 0.0f));
                if(xmlState->hasAttribute("SourceElevDeg" + String(i)))
                    panner_setSourceElev_deg(hPan, i, (float)xmlState->getDoubleAttribute("SourceElevDeg" + String(i), 0.0f));
            }
            if(xmlState->hasAttribute("nSources"))
                panner_setNumSources(hPan, xmlState->getIntAttribute("nSources", 1));
            if(xmlState->hasAttribute("DTT"))
                panner_setDTT(hPan, (float)xmlState->getDoubleAttribute("DTT", 0.5f));
            
            if(xmlState->hasAttribute("JSONFilePath"))
                lastDir = xmlState->getStringAttribute("JSONFilePath", "");
            
            for(int i=0; i<panner_getMaxNumLoudspeakers(); i++){
                if(xmlState->hasAttribute("LoudspeakerAziDeg" + String(i)))
                    panner_setLoudspeakerAzi_deg(hPan, i, (float)xmlState->getDoubleAttribute("LoudspeakerAziDeg" + String(i),0.0f));
                if(xmlState->hasAttribute("LoudspeakerElevDeg" + String(i)))
                    panner_setLoudspeakerElev_deg(hPan, i, (float)xmlState->getDoubleAttribute("LoudspeakerElevDeg" + String(i), 0.0f));
            }
            if(xmlState->hasAttribute("nLoudspeakers"))
                panner_setNumLoudspeakers(hPan, xmlState->getIntAttribute("nLoudspeakers", 1));
            
            panner_refreshSettings(hPan);
        }
    }
}

//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new PluginProcessor();
}

/* Adapted from the AllRADecoder by Daniel Rudrich, (c) 2017 (GPLv3 license) */
void PluginProcessor::saveConfigurationToFile (File destination, int srcOrLs)
{
    DynamicObject* jsonObj;
    char versionString[10];
    elements.removeAllChildren(nullptr);
    switch (srcOrLs){
        case 0:{
            for (int i=0; i<panner_getNumSources(hPan);i++) {
                elements.appendChild (ConfigurationHelper::
                                      createElement(panner_getSourceAzi_deg(hPan, i),
                                                    panner_getSourceElev_deg(hPan, i),
                                                    1.0f, i+1, false, 1.0f), nullptr);
            }
        }
        break;
        case 1:{
            for (int i=0; i<panner_getNumLoudspeakers(hPan);i++) {
                elements.appendChild (ConfigurationHelper::
                                     createElement(panner_getLoudspeakerAzi_deg(hPan, i),
                                                       panner_getLoudspeakerElev_deg(hPan, i),
                                                       1.0f, i+1, false, 1.0f), nullptr);
            }
        }
        break;
    }
    jsonObj = new DynamicObject();
    jsonObj->setProperty("Name", var("SPARTA Panner source/loudspeaker directions."));
    strcpy(versionString, "v");
    strcat(versionString, JucePlugin_VersionString);
    jsonObj->setProperty("Description", var("This configuration file was created with the SPARTA Panner " + String(versionString) + " plug-in. " + Time::getCurrentTime().toString(true, true)));
    jsonObj->setProperty ("GenericLayout", ConfigurationHelper::convertElementsToVar (elements, "Source/Loudspeaker Directions"));
    Result result2 = ConfigurationHelper::writeConfigurationToFile (destination, var (jsonObj));
}

/* Adapted from the AllRADecoder by Daniel Rudrich, (c) 2017 (GPLv3 license) */
void PluginProcessor::loadConfiguration (const File& configFile, int srcOrLs)
{
    int channelIDs[MAX(MAX_NUM_INPUTS, MAX_NUM_OUTPUTS)+1] = {0};
    int virtual_channelIDs[MAX(MAX_NUM_INPUTS, MAX_NUM_OUTPUTS)+1] = {0};
    elements.removeAllChildren(nullptr);
    Result result = ConfigurationHelper::parseFileForGenericLayout(configFile, elements, nullptr);
    if(result.wasOk()){
        int num_el, num_virtual_el, el_idx, jj;
        num_el = num_virtual_el = el_idx = jj = 0;
        /* get Channel IDs and find number of directions and virtual directions */
        for (ValueTree::Iterator it = elements.begin(); it != elements.end(); ++it){
            if ( !((*it).getProperty("Imaginary"))){
                num_el++; channelIDs[jj] = (*it).getProperty("Channel");
            }
            else{
                virtual_channelIDs[num_virtual_el] = (*it).getProperty("Channel");
                num_virtual_el++; channelIDs[jj] = -1;
            }
            jj++;
        }
        /* remove virtual channels and shift the channel indices down */
        for(int i=0; i<num_virtual_el; i++)
            for(int j=0; j<num_el+num_virtual_el; j++)
                if(channelIDs[j] == -1)
                    for(int k=j; k<num_el+num_virtual_el; k++)
                        channelIDs[k] = channelIDs[k+1];
        
        /* then decriment the channel IDs to remove the gaps */
        for(int i=0; i<num_virtual_el; i++)
            for(int j=0; j<num_el+num_virtual_el; j++)
                if( channelIDs[j] > virtual_channelIDs[i]-i )
                    channelIDs[j]--;
        
        /* update with the new configuration  */
        switch(srcOrLs){
            case 0:{
                panner_setNumSources(hPan, num_el);
                for (ValueTree::Iterator it = elements.begin() ; it != elements.end(); ++it){
                    if ( !((*it).getProperty("Imaginary"))){
                        panner_setSourceAzi_deg(hPan, channelIDs[el_idx]-1, (*it).getProperty("Azimuth"));
                        panner_setSourceElev_deg(hPan, channelIDs[el_idx]-1, (*it).getProperty("Elevation"));
                        el_idx++;
                    }
                }
            }
            break;
            case 1:{
                panner_setNumLoudspeakers(hPan, num_el);
                for (ValueTree::Iterator it = elements.begin() ; it != elements.end(); ++it){
                    if ( !((*it).getProperty("Imaginary"))){
                        panner_setLoudspeakerAzi_deg(hPan, channelIDs[el_idx]-1, (*it).getProperty("Azimuth"));
                        panner_setLoudspeakerDist_deg(hPan, channelIDs[el_idx] - 1, (*it).getProperty("Distance")); //XXXX
                        panner_setLoudspeakerElev_deg(hPan, channelIDs[el_idx]-1, (*it).getProperty("Elevation"));
                        el_idx++;
                    }
                }
            }
            break;
        }
    }
}
