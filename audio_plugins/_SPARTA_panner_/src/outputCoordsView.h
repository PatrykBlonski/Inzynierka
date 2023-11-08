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

//[/Headers]



//==============================================================================
/**
                                                                    //[Comments]
    An auto-generated component, created by the Projucer.

    Describe your class and how it works here!
                                                                    //[/Comments]
*/
class outputCoordsView  : public Component
                          
{
public:
    //==============================================================================
    outputCoordsView (PluginProcessor* ownerFilter, int _maxNCH, int _currentNCH );
    ~outputCoordsView() override;

    //==============================================================================
    //[UserMethods]     -- You can add your own custom methods in this section.
    void setNCH(int newNCH){
		newNCH = newNCH > MAX_NUM_OUTPUTS ? MAX_NUM_OUTPUTS : newNCH;
        refreshCoords();
		if (newNCH != currentNCH) {
			currentNCH = newNCH;
			resized();
			labelHasChanged = true;
		}
    }
    bool getHasALabelChanged(){ return labelHasChanged; }
    void setHasALabelChange(bool newState){ labelHasChanged = newState; }

    //[/UserMethods]

    void paint (juce::Graphics& g) override;
    void resized() override;
  //  void labelValueChanged (juce::Label* labelThatWasMoved);



private:
    //[UserVariables]   -- You can add your own custom variables in this section.
    PluginProcessor* hVst;
    void *hPan;
    void refreshCoords();
    std::unique_ptr<Label>* aziLabels;
    std::unique_ptr<Label>* elevLabels;
    std::unique_ptr<Label>* distLabels;
    int maxNCH, currentNCH;
    bool labelHasChanged;
    //[/UserVariables]

    //==============================================================================
    std::unique_ptr<juce::Label> dummyLabel;


    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (outputCoordsView)
};

//[EndFile] You can add extra defines here...
//[/EndFile]

