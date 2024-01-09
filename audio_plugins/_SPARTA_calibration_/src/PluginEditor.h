/*
  ==============================================================================

  This is an automatically generated GUI class created by the Projucer!

  Be careful when adding custom code to these files, as only the code within
  the "//[xyz]" and "//[/xyz]" sections will be retained when the file is loaded
  and re-saved.

  Created with Projucer version: 7.0.7

  ------------------------------------------------------------------------------

  The Projucer is part of the JUCE library.
  Copyright (c) 2020 - Raw Material Software Limited.

  ==============================================================================
*/

#pragma once

//[Headers]     -- You can add your own extra header files here --

#include "JuceHeader.h"
#include "PluginProcessor.h"
#include "outputCoordsView.h"
#include "_common.h"
#include "calibrationView.h"
#include "overlay.h"
#include "../../resources/SPARTALookAndFeel.h"



typedef enum _SPARTA_WARNINGS{
    k_warning_none,
    k_warning_frameSize,
    k_warning_supported_fs,
    k_warning_NinputCH,
    k_warning_NoutputCH
}SPARTA_WARNINGS;
//[/Headers]



//==============================================================================
/**
                                                                    //[Comments]
    An auto-generated component, created by the Introjucer.

    Describe your class and how it works here!
                                                                    //[/Comments]
*/
class PluginEditor  : public AudioProcessorEditor,
                      public MultiTimer,
                      public juce::ComboBox::Listener,
                      public juce::Slider::Listener,
                      public juce::Button::Listener
{
public:
    //==============================================================================
    PluginEditor (PluginProcessor* ownerFilter);
    ~PluginEditor() override;

    //==============================================================================
    //[UserMethods]     -- You can add your own custom methods in this section.


    //[/UserMethods]

    void paint (juce::Graphics& g) override;
    void resized() override;
    void comboBoxChanged (juce::ComboBox* comboBoxThatHasChanged) override;
    void sliderValueChanged (juce::Slider* sliderThatWasMoved) override;
    void buttonClicked (juce::Button* buttonThatWasClicked) override;


private:
    //[UserVariables]   -- You can add your own custom variables in this section.
    PluginProcessor* hVst;
    void* hPan;
    void* hPm;
    int switchOff{};
    void timerCallback(int timerID) override;
#ifndef PLUGIN_EDITOR_DISABLE_OPENGL
    std::unique_ptr<OpenGLGraphicsContextCustomShader> shader;
    OpenGLContext openGLContext;
#endif
    double progress = 0.0;
    ProgressBar progressbar;

    /* Look and Feel */
    SPARTALookAndFeel LAF;

    /* loudspeaker coordinates viewport */
    std::unique_ptr<Viewport> loudspeakerCoordsVP;
    outputCoordsView* loudspeakerCoordsView_handle;

    /* json file loading/saving */
    std::unique_ptr<juce::FileChooser> chooser;

    /* panning window */
    std::unique_ptr<calibrationView> panWindow;
    bool refreshPanViewWindow;

    std::unique_ptr<overlay> overlayWindow;
    bool refreshOverlayWindow;

    /* warnings */
    SPARTA_WARNINGS currentWarning;

    /* tooltips */
    SharedResourcePointer<TooltipWindow> tipWindow;
    std::unique_ptr<juce::ComboBox> pluginDescription; /* Dummy combo box to provide plugin description tooltip */
    std::unique_ptr<juce::ComboBox> windowDescription; /* Dummy combo box to provide plugin description tooltip */

    //[/UserVariables]

    //==============================================================================
   // std::unique_ptr<juce::ComboBox> CBsourceDirsPreset;
    std::unique_ptr<juce::Slider> SL_num_sources;
    std::unique_ptr<juce::ToggleButton> TB_showOutputs;
    std::unique_ptr<juce::Slider> SL_num_loudspeakers;
    //std::unique_ptr<juce::TextButton> tb_loadJSON_src;
    std::unique_ptr<juce::TextButton> tb_saveJSON_ls;
    std::unique_ptr<juce::TextButton> tb_calibration;
    std::unique_ptr<juce::ComboBox> CBformat;
    std::unique_ptr<juce::ComboBox> CBnorm;
    std::unique_ptr<juce::TextButton> tb_play;
    std::unique_ptr<juce::TextButton> tb_next;
    std::unique_ptr<juce::ToggleButton> TB_PlayAll;



    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (PluginEditor)

};

//[EndFile] You can add extra defines here...
//[/EndFile]

