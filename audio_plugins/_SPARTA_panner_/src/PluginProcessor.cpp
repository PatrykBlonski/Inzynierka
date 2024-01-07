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
    panner_create(&hPan);
    powermap_create(&hPm);
    TempBuffer.setSize(4, 1024);
    refreshWindow = true;
    startTimer(TIMER_PROCESSING_RELATED, 80); 
}

PluginProcessor::~PluginProcessor()
{
    panner_destroy(&hPan);
    powermap_destroy(&hPm);
}

void PluginProcessor::generateSineSweep(float sampleRate, juce::AudioBuffer<float>& sweepBuffer) {
    // Calculate the total number of samples
    int totalSamples = static_cast<int>(duration * sampleRate);

    // Resize the sweepBuffer to hold the sine sweep
    sweepBuffer.setSize(1, totalSamples);  // Assuming mono

    // Variables for sweep generation
    double phase = 0.0;
    double phaseIncrement;
    double frequency;
    double time = 0.0;
    double timeIncrement = 1.0 / sampleRate;

    for (int i = 0; i < totalSamples; ++i) {
        // Compute the instantaneous frequency using logarithmic interpolation
        frequency = startFrequency * std::pow(endFrequency / startFrequency, time / duration);

        // Compute the phase increment
        phaseIncrement = 2.0 * double_Pi * frequency / sampleRate;

        // Generate the sine sample
        sweepBuffer.setSample(0, i, std::sin(phase));

        // Update phase and wrap around 2*PI if necessary
        phase += phaseIncrement;
        if (phase >= 2.0 * double_Pi) {
            phase -= 2.0 * double_Pi;
        }

        // Increment time
        time += timeIncrement;
    }
}


void PluginProcessor::startCalibration() {
    generateSineSweep(48000, SweepBuffer);
   // juce::String juceStringPath = "D:\\STUDIA\\7sem\\impulse_responses\\source.wav";
   //// std::string stdStringPath = "C:\\Users\\patry\\Downloads\\sweep.wav";
   // SweepBuffer.setSize(1, 480000+256);
   // SweepBuffer = loadImpulseResponse(juceStringPath);
    dists.clear();
    loudspeakerNumber = 0;
    juce::String fileName = "recorded_audio_sweep_" + juce::String(loudspeakerNumber) + ".wav";
    saveBufferToWav(SweepBuffer, fileName);
    prepareToPlay(getSampleRate(), 480000);
    recordingBuffer.setSize(panner_getNumSources(hPan), getSampleRate() * duration);
    recordingBuffer.clear();
    latency = 0;
     juce::String juceStringPath = "D:\\STUDIA\\7sem\\impulse_responses\\conv_signal_" + juce::String(loudspeakerNumber) + ".wav";
    //juceStringPath = "D:\\STUDIA\\7sem\\impulse_responses\\conv_signal_" + juce::String(loudspeakerNumber) + ".wav";
    ImpulseBuffer = loadImpulseResponse(juceStringPath);
    currentRecordingPosition = 0;
    calibrating = true;
    //phase = 0.0;
    //timeElapsed = 0.0;
    //frequency = startFrequency;
    isRecording = true;
}

void PluginProcessor::endCalibration() {
    isRecording = false;
    calibrating = false;
   
    juce::String fileName = "recorded_audio_mic_" + juce::String(loudspeakerNumber) + ".wav";
    saveBufferToWav(recordingBuffer, fileName);
   /* fileName = "recorded_audio_sweep_" + juce::String(loudspeakerNumber) + ".wav";
    saveBufferToWav(SweepBuffer, fileName);*/
    distanceCalculation(SweepBuffer, ImpulseBuffer, loudspeakerNumber);
    panner_process(hPan, ImpulseBuffer.getArrayOfWritePointers(), panner_getNumSources(hPan), ImpulseBuffer.getNumSamples(), loudspeakerNumber, 1);
    refreshWindow = true;

    loudspeakerNumber++;
    if (loudspeakerNumber < panner_getNumLoudspeakers(hPan)) {
        latency = 0;
        if (playAll) {
            ImpulseBuffer.clear();
            juce::String juceStringPath = "D:\\STUDIA\\7sem\\impulse_responses\\conv_signal_" + juce::String(loudspeakerNumber) + ".wav";
            ImpulseBuffer = loadImpulseResponse(juceStringPath);
            currentRecordingPosition = 0;
            calibrating = true;
        }
       
        //phase = 0.0;
        //timeElapsed = 0.0;
        //frequency = startFrequency;
    }
}

void PluginProcessor::next() {
    currentRecordingPosition = 0;
    latency = 0;

    ImpulseBuffer.clear();
    juce::String juceStringPath = "D:\\STUDIA\\7sem\\impulse_responses\\conv_signal_" + juce::String(loudspeakerNumber) + ".wav";
    ImpulseBuffer = loadImpulseResponse(juceStringPath);
    recordingBuffer.clear();

    calibrating = true;

}

double PluginProcessor::computeSweepFrequency(double time) {
    return jmap(time, 0.0, duration, startFrequency, endFrequency);        
}

void PluginProcessor::setParameter (int index, float newValue)
{
    /* standard parameters */
    if(index < k_NumOfParameters){
        switch (index) {
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
            newValueScaled = newValue * 10.0f;
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
    powermap_init(hPm, nSampleRate);
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
        //frequency = computeSweepFrequency(timeElapsed);

        //// Compute the sample of the sine sweep
        //float sineSample = std::sin(phase);

        //// Write this sample to all output channels (you can change this if needed)
        //buffer.setSample(loudspeakerNumber + numberOfInputs, sample, sineSample);

        //// Increment phase and wrap if it exceeds 2*PI
        //phase += (2.0 * double_Pi * frequency) * deltaT;
        //if (phase > 2.0 * double_Pi) {
        //    phase -= 2.0 * double_Pi;
        //}

        //// Update time elapsed
        //timeElapsed += deltaT;

        //// Stop the sweep if the duration has been reached
        //if (timeElapsed >= duration) {
        //    endCalibration();
        //    break;
        //}
        buffer.setSample(numberOfInputs + loudspeakerNumber, sample, SweepBuffer.getSample(0, currentRecordingPosition + sample));
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
    float** bufferData;


    // Clear any unused channels
    for (int i = totalNumInputChannels; i < totalNumOutputChannels; ++i) {
        buffer.clear(i, 0, numSamples);
    }


    if (calibrating) {
            if (currentRecordingPosition + buffer.getNumSamples() <= recordingBuffer.getNumSamples()) 
                buffer.copyFrom(numberOfInputs + loudspeakerNumber, 0, SweepBuffer, 0, currentRecordingPosition, buffer.getNumSamples());
            if (currentRecordingPosition - 1 * buffer.getNumSamples() <= recordingBuffer.getNumSamples() && latency > 1) {
                for (int channel = 0; channel < numberOfInputs; ++channel) {
                    recordingBuffer.copyFrom(channel, currentRecordingPosition - 2 * buffer.getNumSamples(), buffer, channel, 0, buffer.getNumSamples());
                }
            }
            else if (currentRecordingPosition - 1 * buffer.getNumSamples() > recordingBuffer.getNumSamples()){
                endCalibration();
            }
            currentRecordingPosition += buffer.getNumSamples();
            latency += 1;
    }
    if (isPlaying) {
        if (currentRecordingPosition + 1024 <= ImpulseBuffer.getNumSamples()) {
            for (int i = 0; i < 4; i++) {
                TempBuffer.copyFrom(i, 0, ImpulseBuffer, i, currentRecordingPosition, 1024);
            }
            bufferData = TempBuffer.getArrayOfWritePointers();
            currentRecordingPosition+=1024;
            powermap_analysis(hPm, bufferData, 4, 1024, isPlaying);
        }
        else {
            isPlaying = false;
        }
    }
}

void PluginProcessor::distanceCalculation(AudioBuffer<float>& sweep, AudioBuffer<float>& input, int loudNum)
{
    const int numOfSamples = sweep.getNumSamples();
    const int fftSize = nextPowerOfTwo(numOfSamples); //finding next power of 2 for FFT calc

    kiss_fft_cfg fft_cfg = kiss_fft_alloc(fftSize, 0, nullptr, nullptr);
    kiss_fft_cfg ifft_cfg = kiss_fft_alloc(fftSize, 1, nullptr, nullptr);

    std::vector<kiss_fft_cpx> kiss_sweep_in(fftSize);
    std::vector<kiss_fft_cpx> kiss_sweep_out(fftSize);
    std::vector<kiss_fft_cpx> kiss_input_in(fftSize);
    std::vector<kiss_fft_cpx> kiss_input_out(fftSize);

    for (int i = 0; i < numOfSamples; ++i) {
        kiss_sweep_in[i].r = sweep.getSample(0, i);
        kiss_sweep_in[i].i = 0.0f;

        kiss_input_in[i].r = input.getSample(0, i);
        kiss_input_in[i].i = 0.0f;
    }
    // Zero-pad the rest if needed
    for (int i = numOfSamples; i < fftSize; ++i) {
        kiss_sweep_in[i].r = 0.0f;
        kiss_sweep_in[i].i = 0.0f;

        kiss_input_in[i].r = 0.0f;
        kiss_input_in[i].i = 0.0f;
    }

    kiss_fft(fft_cfg, kiss_sweep_in.data(), kiss_sweep_out.data());
    kiss_fft(fft_cfg, kiss_input_in.data(), kiss_input_out.data());

    const float epsilon = 1e-3f; // Threshold to prevent division by very small numbers

    for (int i = 0; i < fftSize; ++i) {
        if (std::abs(kiss_sweep_out[i].r) > epsilon || std::abs(kiss_sweep_out[i].i) > epsilon) {
            // Perform the division
            kiss_input_in[i].r = (kiss_input_out[i].r * kiss_sweep_out[i].r + kiss_input_out[i].i * kiss_sweep_out[i].i) /
                (kiss_sweep_out[i].r * kiss_sweep_out[i].r + kiss_sweep_out[i].i * kiss_sweep_out[i].i);
            kiss_input_in[i].i = (kiss_input_out[i].i * kiss_sweep_out[i].r - kiss_input_out[i].r * kiss_sweep_out[i].i) /
                (kiss_sweep_out[i].r * kiss_sweep_out[i].r + kiss_sweep_out[i].i * kiss_sweep_out[i].i);
        }
        else {
            // Regularization
            kiss_input_in[i].r = 0.0f;
            kiss_input_in[i].i = 0.0f;
        }
    }

    kiss_fft(ifft_cfg, kiss_input_in.data(), kiss_input_out.data());

    // Create an AudioBuffer to store the impulse response
    AudioBuffer<float> impulseResponse(1, fftSize);

    //Variables for fiding time
    float peakTime = 0.0f;
    float maxVal = std::numeric_limits<float>::lowest();
    int peakIndex = 0;  // Index of the peak value

    for (int i = 0; i < fftSize; ++i) {
        impulseResponse.setSample(0, i, kiss_input_out[i].r);

        float currentVal = kiss_input_out[i].r;
        if (std::abs(currentVal) > maxVal)
        {
            maxVal = std::abs(currentVal);
            peakIndex = i; // Update index of the peak value
        }
    }

    kiss_fft_free(fft_cfg);
    kiss_fft_free(ifft_cfg);

    juce::String fileName = "impulse_" + juce::String(loudspeakerNumber) + ".wav";
    saveBufferToWav(impulseResponse, fileName);

    peakTime = static_cast<float>(peakIndex) / static_cast<float>(48000); // Use actual sample rate

    const float speedOfSound = 343.0f;
    float distance = peakTime * speedOfSound;

    //panner_setLoudspeakerDist_deg(hPan, loudNum, distance);
    dists.push_back(distance);
    float avgDist = 0;
    int i{};
    std::vector<float> norm;
    float max = *std::max_element(dists.begin(), dists.end());
    for (auto& elem : dists) {
        avgDist += elem / (loudspeakerNumber + 1);
        norm.push_back(sqrt(elem / max));
    }
    float maxNorm = *std::max_element(norm.begin(), norm.end());
    float minNorm = *std::min_element(norm.begin(), norm.end());
    for (auto& elem : norm) {
        panner_setLoudspeakerDist_plot(hPan, i, 0.7 + 0.3 * (elem - minNorm) / (maxNorm - minNorm));
        panner_setLoudspeakerDist_deg(hPan, i, dists[i] - avgDist);
        i++;
    }
    
   // // Check if buffers are valid and have the same number of samples
   // jassert(sweep.getNumSamples() == input.getNumSamples());
   //// jassert(sweep.getNumChannels() > 0 && input.getNumChannels() > 0);

   // const int numOfSamples = sweep.getNumSamples();
   // const int fftSize = nextPowerOfTwo(numOfSamples); //finding next power of 2 for FFT calc
   // juce::dsp::FFT fft((int)std::log2(fftSize));         

   // // Allocate space for the FFT data
   // HeapBlock<float> sweepFFTData(fftSize * 2, true); // Initialize to zero
   // HeapBlock<float> inputFFTData(fftSize * 2, true); // Initialize to zero

   // // Copy input data to the fftData's real part
   // const float* sweepSamples = sweep.getReadPointer(0);
   // const float* inputSamples = input.getReadPointer(0);

   // // Prepare the data for FFT
   // for (int i = 0; i < numOfSamples; ++i)
   // {
   //     sweepFFTData[i] = sweepSamples[i];
   //     inputFFTData[i] = inputSamples[i];
   // }

   // // Perform the FFT in-place
   // fft.performFrequencyOnlyForwardTransform(sweepFFTData);
   // fft.performFrequencyOnlyForwardTransform(inputFFTData);

   // // Allocate space for the result of the division (transfer function)
   // HeapBlock<juce::dsp::Complex<float>> transferFunctionFFTData(fftSize, true); // Initialize to zero
   // 
   // // Perform the division in the frequency domain with regularization
   // const float epsilon = 1e-3f; // Threshold to prevent division by very small numbers
   // for (int i = 0; i < fftSize; ++i)
   // {
   //     juce::dsp::Complex<float> sweepValue(sweepFFTData[i * 2], sweepFFTData[i * 2 + 1]);
   //     juce::dsp::Complex<float> inputValue(inputFFTData[i * 2], inputFFTData[i * 2 + 1]);

   //     if (std::abs(sweepValue) > epsilon)
   //     {
   //         transferFunctionFFTData[i] = inputValue / sweepValue;
   //     }
   //     else
   //     {
   //         transferFunctionFFTData[i] = juce::dsp::Complex<float>(0.0f, 0.0f); // Regularization
   //     }
   // }

   // // Prepare buffer for inverse FFT result (time domain)
   // HeapBlock<float> estimatedIRData(fftSize * 2, true); // Initialize to zero

   // // Prepare data for inverse FFT
   // for (int i = 0; i < fftSize; ++i)
   // {
   //     estimatedIRData[i * 2] = transferFunctionFFTData[i].real();
   //     estimatedIRData[i * 2 + 1] = transferFunctionFFTData[i].imag();
   // }

   // // Perform the inverse FFT
   // fft.performRealOnlyInverseTransform(estimatedIRData);

   // // Find the peak in the impulse response to determine the time of arrival
   // float peakTime = 0.0f;
   // float maxVal = std::numeric_limits<float>::lowest();
   // AudioBuffer<float> impulseResponse;
   // impulseResponse.setSize(1, fftSize);

   // for (int i = 1; i < fftSize - 1; ++i)
   // {
   //     float val = estimatedIRData[i * 2]; // We only need the real part
   //     impulseResponse.setSample(0, i - 1, val);
   //     if (val > maxVal)
   //     {
   //         maxVal = val; //finding peak Value
   //         peakTime = (float)i / 48000.0f; //const sampling rate
   //     }
   // }
   // juce::String fileName = "impulse_" + juce::String(loudspeakerNumber) + ".wav";
   // saveBufferToWav(impulseResponse, fileName);
   // // Calculate the distance using the time of arrival and the speed of sound
   // const float speedOfSound = 343.0f; // Speed of sound in m/s at 20 degrees Celsius
   // float distance = peakTime * speedOfSound; //results 2.2m-2.6m
   // panner_setLoudspeakerDist_deg(hPan, loudNum, distance);
}


void PluginProcessor::loadImpulseResponse()
{
    // Assuming you have a method to load an AudioBuffer with your IR...
    std::string stdStringPath = "C:\\Users\\patry\\Documents\\impulse_response_5.wav";
    //std::string stdStringPath = "C:\\Users\\patry\\Downloads\\sweep.wav";
    juce::String juceStringPath(stdStringPath);
    AudioBuffer<float> irBuffer = loadImpulseResponse(juceStringPath);
   // saveBufferToWav(irBuffer);
    // Load the impulse response for each channel into each convolution instance
    if (irBuffer.getNumChannels() == 4)
    {
        // You should create an AudioBuffer for each channel
        AudioBuffer<float> singleChannelIR1(1, irBuffer.getNumSamples());
        AudioBuffer<float> singleChannelIR2(1, irBuffer.getNumSamples());
        AudioBuffer<float> singleChannelIR3(1, irBuffer.getNumSamples());
        AudioBuffer<float> singleChannelIR4(1, irBuffer.getNumSamples());
        singleChannelIR1.copyFrom(0, 0, irBuffer, 0, 0, irBuffer.getNumSamples());
        singleChannelIR2.copyFrom(0, 0, irBuffer, 1, 0, irBuffer.getNumSamples());
        singleChannelIR3.copyFrom(0, 0, irBuffer, 2, 0, irBuffer.getNumSamples());
        singleChannelIR4.copyFrom(0, 0, irBuffer, 3, 0, irBuffer.getNumSamples());
        // Load the impulse response into the convolution object
        convolution1.loadImpulseResponse(
            std::move(singleChannelIR1),
            getSampleRate(),
            dsp::Convolution::Stereo::no,
            dsp::Convolution::Trim::no,
            dsp::Convolution::Normalise::yes
        );

        convolution2.loadImpulseResponse(
            std::move(singleChannelIR2),
            getSampleRate(),
            dsp::Convolution::Stereo::no,
            dsp::Convolution::Trim::no,
            dsp::Convolution::Normalise::yes
        );

        convolution3.loadImpulseResponse(
            std::move(singleChannelIR3),
            getSampleRate(),
            dsp::Convolution::Stereo::no,
            dsp::Convolution::Trim::no,
            dsp::Convolution::Normalise::yes
        );

        convolution4.loadImpulseResponse(
            std::move(singleChannelIR4),
            getSampleRate(),
            dsp::Convolution::Stereo::no,
            dsp::Convolution::Trim::no,
            dsp::Convolution::Normalise::yes
        );
        // ... Repeat for convolution2, convolution3, convolution4 with appropriate channels
    }
    else
    {
        DBG("The impulse response buffer does not have enough channels.");
    }
}

void PluginProcessor::saveBufferToWav(AudioBuffer<float> &buffer, const String& path)
{
    // File path

    juce::File outputFile = juce::File::getSpecialLocation(juce::File::userDesktopDirectory).getChildFile(path);

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
        buffer.getNumChannels(),
        32,  // Bit depth, e.g., 16 bits
        {},  // No metadata
        0    // Use the default quality
    ));

    if (writer.get() != nullptr)
    {
        // Write the entire buffer to the file
        writer->writeFromAudioSampleBuffer(buffer, 0, buffer.getNumSamples());
    }
}

AudioBuffer<float> PluginProcessor::loadImpulseResponse(const String& filePath)
{
    AudioBuffer<float> buffer;

    File file(filePath);
    if (!file.existsAsFile())
    {
        DBG("File does not exist: " << filePath);
        return buffer;
    }

    AudioFormatManager formatManager;
    formatManager.registerBasicFormats();
    auto variable = formatManager.createReaderFor(file);

    std::unique_ptr<AudioFormatReader> reader(variable);
    if (reader)
    {
        if (reader->usesFloatingPointData) // Check if the reader can handle floating-point data
        {
            buffer.setSize((int)reader->numChannels, (int)reader->lengthInSamples);
            reader->read(&buffer, 0, (int)reader->lengthInSamples, 0, true, true);
        }
        else
        {
            DBG("The file is not in a floating-point format.");
        }
    }
    else
    {
        DBG("Failed to create audio format reader for file: " << filePath);
    }

    return buffer;
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
                                                       panner_getLoudspeakerDist_deg(hPan, i), i + 1, false, 1.0f), nullptr);
            }
        }
        break;
    }
    jsonObj = new DynamicObject();
    jsonObj->setProperty("Name", var("SPARTA Calibration loudspeaker directions."));
    strcpy(versionString, "v");
    strcat(versionString, JucePlugin_VersionString);
    jsonObj->setProperty("Description", var("This configuration file was created with the SPARTA Calibration " + String(versionString) + " plug-in. " + Time::getCurrentTime().toString(true, true)));
    jsonObj->setProperty ("GenericLayout", ConfigurationHelper::convertElementsToVar (elements, "Loudspeaker Directions"));
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
