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
class LoudspeakerIconComponent : public Component,
    public SettableTooltipClient {
public:
    //LoudspeakerIconComponent(float distance, float azimuth, float elevation)
    //    : distanceValue(distance), azimuthValue(azimuth), elevationValue(elevation) {
    //    // Set the tooltip text for this component
    //}

    void paint(Graphics& g) override {
        g.setColour(Colour::fromFloatRGBA(0.5f, 1.0f, 0.1f, 1.0f));
        g.setOpacity(0.3f);
        g.fillRect(getLocalBounds()); // Draw the rectangle for the icon

    }

    // You can override mouseEnter and mouseExit if you want to change the appearance on hover
    // void mouseEnter(const MouseEvent& event) override { ... }
    // void mouseExit(const MouseEvent& event) override { ... }

    // Add setters if you need to change the values dynamically
    void setDistance(float newDistance) {
        distanceValue = newDistance;
        // Update the tooltip or other properties if needed
    }
    void setTooltip(const String& newTooltip) override;
    // Add other setters for azimuth and elevation if needed

private:
    float distanceValue;
    float azimuthValue;
    float elevationValue;

    static constexpr int icon_size = 10; // Example icon size, adjust as needed
};


//==============================================================================
/**
                                                                    //[Comments]
    An auto-generated component, created by the Projucer.

    Describe your class and how it works here!
                                                                    //[/Comments]
*/
class calibrationView  : public Component
{
public:
    //==============================================================================
    calibrationView (PluginProcessor* ownerFilter, int _width, int _height);
    ~calibrationView() override;

    //==============================================================================
    //[UserMethods]     -- You can add your own custom methods in this section.

    void refreshPanView();
    void setShowInputs(bool state){ showInputs = state; }
    void setShowOutputs(bool state){ showOutputs = state; }
    bool getSourceIconIsClicked(){ return sourceIconIsClicked; }

    //[/UserMethods]

    void paint (juce::Graphics& g) override;
    void resized() override;
    void mouseDown (const juce::MouseEvent& e) override;
    void mouseDrag (const juce::MouseEvent& e) override;
    void mouseUp (const juce::MouseEvent& e) override;



private:
    //[UserVariables]   -- You can add your own custom variables in this section.
    PluginProcessor* hVst;
    void* hPan;
    int width;
    int height;
    bool showInputs;
    bool showOutputs;
    Rectangle<float> SourceIcons[MAX_NUM_INPUTS];
    LoudspeakerIconComponent* LoudspeakerIcons[MAX_NUM_OUTPUTS];
    //Rectangle<float> LoudspeakerIcons[MAX_NUM_OUTPUTS];
    int NSources;
    int NLoudspeakers;
    int NLoudspeakersPrev;
    bool sourceIconIsClicked;
    int indexOfClickedSource;
    //[/UserVariables]

    //==============================================================================


    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (calibrationView)
};

//[EndFile] You can add extra defines here...
//[/EndFile]

