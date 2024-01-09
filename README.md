# SPARTA

Spatial Audio Real-Time Applications (SPARTA) [1]. A collection of VST audio plug-ins for spatial audio production, reproduction and visualisation. Developed using [JUCE](https://github.com/WeAreROLI/JUCE/) and the [Spatial_Audio_Framework](https://github.com/leomccormack/Spatial_Audio_Framework).

## Plug-in descriptions
* **Calibration** - A calibrating plugin, estimating location of sound sources

## Pre-built plug-ins

audio_plugins/_SPARTA_calibration_/Builds/VisualStudio2022/x64/Debug/VST/sparta_calibration.dll

## Building the plug-ins yourself

First clone the repository (including submodules) with:

```
git clone --recursive https://github.com/leomccormack/SPARTA
# or if you have already cloned the repository, update/init with:
git submodule update --init --recursive
```

## Prerequisites 

The [VST2_SDK](https://web.archive.org/web/20181016150224/https://download.steinberg.net/sdk_downloads/vstsdk3610_11_06_2018_build_37.zip) should be placed in the 'SDKs' folder like so:
```
SDKs/VST2_SDK
```

By default, **MacOSX, Linux and Windows (x86_64/amd64)** users need to install [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html) (MKL and IPP) and run the **install-safmkl**.sh/.bat and **install-safipp**.sh/.bat scripts found in SDKs/Spatial_Audio_Framework/scripts. Whereas, **Raspberry Pi (ARM)** users instead require OpenBLAS and LAPACKE libraries:
``` 
sudo apt-get install liblapack3 liblapack-dev libopenblas-base libopenblas-dev liblapacke-dev
```
Note, however, that alternative performance libraries may also be used, with more information provided [here](https://github.com/leomccormack/Spatial_Audio_Framework/blob/master/docs/PERFORMANCE_LIBRARY_INSTRUCTIONS.md).

**Linux (x86_64/amd64 and ARM)** users must also install the following libraries required by JUCE:

```
sudo apt-get install x11proto-xinerama-dev libwebkit2gtk-4.0-dev libgtk-3-dev x11proto-xext-dev libcurl4-openssl-dev libasound2-dev
```

## Building the plug-ins via CMake 

The plug-ins may be built with CMake (version 3.15 or higher):
 ```
 mkdir build
 cmake -S . -B build -DSAF_ENABLE_SOFA_READER_MODULE=1
 cd build
 make
 ```
 
Or for Visual Studio users (using x64 Native Tools Command Prompt as **administrator**):
```
cmake -S . -B build -G "Visual Studio 15 Win64" -DSAF_ENABLE_SOFA_READER_MODULE=1 
cd build
msbuild ALL_BUILD.vcxproj /p:Configuration=Release /m
```

## Building the plug-ins via the included scripts

**MacOSX/Linux users** may run the following bash script via the Terminal to build all of the plugins:

```
./build-plugins.sh all
# Note: MacOSX users may need to first install and enable Xcode Command Line Tools:
xcode-select --install 
sudo xcode-select -s /Applications/Xcode.app/Contents/Developer 
```

**Windows users** may instead run the following batch script via the "x64 Developer Command Prompt for VS.exe":

```
build-plugins.bat <path/to/Projucer.exe>
```

### Additional scripts and options for MacOSX/Linux users

The repository also includes the following install scripts:
```
./install-juce.sh      # builds a GPLv3 version of the Projucer App and copies it into "SDKs"
./install-vst2_sdk.sh  # downloads, unzips, and places the VST2_SDK into "SDKs"
```

The build.plugins.sh script also supports many additional options:
```
./build-plugins.sh --help    # help information
./build-plugins.sh projuce   # generates Linux makefiles and IDE project files for all plugins
./build-plugins.sh clean     # cleans all plugins 
./build-plugins.sh build     # builds all plugins
./build-plugins.sh all       # projuces, cleans, and then builds all plugins
./build-plugins.sh _SPARTA_ambiBIN_ all       # projuces+cleans+builds sparta_ambiBIN.vst
./build-plugins.sh _SPARTA_ambiENC_ build     # builds "sparta_ambiENC.vst"
./build-plugins.sh _SPARTA_array2sh_ projucer # opens "sparta_array2sh.jucer" with Projucer
```

## Building the plug-ins without scripts or CMake

You may also manually open each .jucer file with the Projucer App and click "Save Project". This will generate Visual Studio (2015/2017) solution files, Xcode project files, Linux Makefiles (amd64), and Raspberry Pi Linux Makefiles (ARM), which are placed in:

```
audio_plugins/_SPARTA_X_/make/
```

To generate project files for other IDEs, you may open and configure the included .jucer files accordingly.

