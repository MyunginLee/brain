// MUS109IA & MAT276IA.
// Spring 2022
// Course Instrument 06. FM Vib-Visual (Object & Spectrum)
// Press '[' or ']' to turn on & off GUI
// Able to play with MIDI device
// ***MacOS may require manual installation of assets3D.
// brew install assimp
// https://github.com/assimp/assimp/blob/master/Build.md
// Myungin Lee

#include <cstdio> // for printing to stdout

#include "Gamma/Analysis.h"
#include "Gamma/Effects.h"
#include "Gamma/Envelope.h"
#include "Gamma/Gamma.h"
#include "Gamma/Oscillator.h"
#include "Gamma/Types.h"
#include "Gamma/DFT.h"

#include "al/app/al_App.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "al/scene/al_PolySynth.hpp"
#include "al/scene/al_SynthSequencer.hpp"
#include "al/ui/al_ControlGUI.hpp"
#include "al/ui/al_Parameter.hpp"
#include "al/math/al_Random.hpp"
#include "al/sound/al_Reverb.hpp"
#include "al_ext/assets3d/al_Asset.hpp"
#include <algorithm>
#include <cstdint>
#include <vector>
#include "soundSynth.cpp"

using namespace gam;
using namespace al;
using namespace std;
#define FFT_SIZE 4048

class MyApp : public App
{
public:
  SynthGUIManager<FMWT> synthManager{"synth6"};
  int midiNote;
  FMWT fm;
  Mesh mSpectrogram;
  vector<float> spectrum;
  bool showGUI = true;
  bool half_draw = false;
  bool showSpectro = true;
  bool navi = false;
  gam::STFT stft = gam::STFT(FFT_SIZE, FFT_SIZE / 4, 0, gam::HANN, gam::MAG_FREQ);
  Scene *ascene1{nullptr};
  Scene *ascene2{nullptr};
  Vec3f scene_min, scene_max, scene_center;
  Vec3f scene_min_half, scene_max_half, scene_center_half;
  vector<Mesh> mMesh;
  vector<Mesh> mMesh_half;
  Texture transparent_tex;

  vector<Mesh> importOBJ(std::string fileName, vector<Mesh> mMesh){
    Scene *ascene{nullptr};
    ascene = Scene::import(fileName);
    if (!ascene)
    {
      printf("error reading %s\n", fileName.c_str());
      // return;
    }
    ascene->getBounds(scene_min, scene_max);
    scene_center = (scene_min + scene_max) / 2.f;
    // ascene->print();
    mMesh.resize(ascene->meshes());
    for (int i = 0; i < ascene->meshes(); i += 1)
    {
      ascene->mesh(i, mMesh[i]);
    }
    return mMesh;
  }

  virtual void onInit() override
  {
    nav().pos(0, 30, 150);
    lens().near(0.1).far(1000).fovy(90);

    // Generate a texture with an alpha channel for transparency
    transparent_tex.create2D(256, 256, Texture::RGBA8);
    int Nx = transparent_tex.width();
    int Ny = transparent_tex.height();
    vector<Colori> pix;
    pix.resize(Nx * Ny);
    std::string fileName1 = "../obj/brain.obj";
    std::string fileName2 = "../obj/brain_half.obj";

    importOBJ(fileName1, mMesh);
    importOBJ(fileName2, mMesh_half);
    
    // Declare the size of the spectrum
    spectrum.resize(FFT_SIZE / 2 + 1);
    imguiInit();
    navControl().active(false); // Disable navigation via keyboard, since we
                                // will be using keyboard for note triggering
    // Set sampling rate for Gamma objects from app's audio
    gam::sampleRate(audioIO().framesPerSecond());
  }

  void onCreate() override
  {
    navControl().active(false);
    // Play example sequence. Comment this line to start from scratch
    //    synthManager.synthSequencer().playSequence("synth2.synthSequence");
    synthManager.synthRecorder().verbose(true);
  }

  void onSound(AudioIOData &io) override
  {
    synthManager.render(io); // Render audio
    while (io())
    {
      io.out(0) = tanh(io.out(0));
      io.out(1) = tanh(io.out(1));
      if (stft(io.out(0)))
      { // Loop through all the frequency bins
        for (unsigned k = 0; k < stft.numBins(); ++k)
        {
          // Here we simply scale the complex sample
          spectrum[k] = tanh(pow(stft.bin(k).real(), 1.3));
          // spectrum[k] = stft.bin(k).real();
        }
      }
    }
  }

  void onAnimate(double dt) override
  {
    navControl().active(navi); // Disable navigation via keyboard, since we
    // Draw GUI
    imguiBeginFrame();
    synthManager.drawSynthControlPanel();
    imguiEndFrame();
    // Map table number to table in memory
    fm.mtable = int(synthManager.voice()->getInternalParameterValue("table"));
  }

  void onDraw(Graphics &g) override
  {
    g.clear(0);
    synthManager.render(g);
    // // Draw Spectrum
    mSpectrogram.reset();
    mSpectrogram.primitive(Mesh::LINE_STRIP);
    if (showSpectro)
    {
      for (int i = 0; i < FFT_SIZE / 2; i++)
      {
        mSpectrogram.color(HSV(0.5 - spectrum[i] * 100));
        mSpectrogram.vertex(i, spectrum[i], 0.0);
      }
      g.meshColor(); // Use the color in the mesh
      g.pushMatrix();
      g.translate(-3, -3, -17);
      g.scale(20.0 / FFT_SIZE, 100, 1.0);
      g.draw(mSpectrogram);
      g.popMatrix();
    }
    gl::depthTesting(true);
    g.lighting(true);
    float tmp = scene_max[0] - scene_min[0];
    // // half brain draw
    if (half_draw)
    {
      g.pushMatrix();
      tmp = scene_max_half[0] - scene_max_half[0];
      tmp = std::max(scene_max_half[1] - scene_max_half[1], tmp);
      tmp = std::max(scene_max_half[2] - scene_max_half[2], tmp);
      tmp = 2.f / tmp;
      g.scale(tmp);
      g.scale(100);
      g.color(RGB(242 / 255., 174 / 255., 177 / 255.));
      for (auto &m : mMesh_half)
      {
        g.draw(m);
      }
      // g.blending(false);
      // g.blendTrans();

      g.popMatrix();
    }
    else    // Full brain draw
    {
      g.pushMatrix();
      tmp = std::max(scene_max[1] - scene_min[1], tmp);
      tmp = std::max(scene_max[2] - scene_min[2], tmp);
      tmp = 2.f / tmp;
      g.scale(tmp);
      g.scale(100);
      g.meshColor();
      g.color(RGB(242 / 255., 174 / 255., 177 / 255.));
      // g.texture();
      // transparent_tex.bind();
      for (auto &m : mMesh)
      {
        g.draw(m);
      }
      g.blending(true);
      g.blendTrans();
      // g.quad(transparent_tex, -1, -1, 2, 2);
      // transparent_tex.unbind();
      g.blending(false);
      g.popMatrix();
    }
    // GUI is drawn here
    if (showGUI)
    {
      imguiDraw();
    }
  }

  bool onKeyDown(Keyboard const &k) override
  {
    if (ParameterGUI::usingKeyboard())
    { // Ignore keys if GUI is using them
      return true;
    }
    if (!navi)
    {
      if (k.shift())
      {
        // If shift pressed then keyboard sets preset
        int presetNumber = asciiToIndex(k.key());
        synthManager.recallPreset(presetNumber);
      }
      else
      {
        // Otherwise trigger note for polyphonic synth
        int midiNote = asciiToMIDI(k.key());
        if (midiNote > 0)
        {
          synthManager.voice()->setInternalParameterValue(
              "frequency", ::pow(2.f, (midiNote - 69.f) / 12.f) * 432.f);
          synthManager.voice()->setInternalParameterValue("table", fm.mtable);
          synthManager.triggerOn(midiNote);
        }
      }
    }
    switch (k.key())
    {
    case ']':
      showGUI = !showGUI;
      break;
    case '[':
      showSpectro = !showSpectro;
      break;
    case '=':
      navi = !navi;
      break;
    case '-':
      half_draw = !half_draw;
      break;
    }
    return true;
  }

  bool onKeyUp(Keyboard const &k) override
  {
    int midiNote = asciiToMIDI(k.key());
    if (midiNote > 0)
    {
      synthManager.triggerOff(midiNote);
    }
    return true;
  }

  void onExit() override { imguiShutdown(); }
};

int main()
{
  MyApp app;

  // Set up audio
  app.configureAudio(48000., 512, 2, 0);

  app.start();
}
