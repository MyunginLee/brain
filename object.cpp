// Brain test
// proof of concept 2023
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
// #include "al/ui/al_ControlGUI.hpp"
#include "al/ui/al_Parameter.hpp"
#include "al/math/al_Random.hpp"
#include "al/sound/al_Reverb.hpp"
#include "al_ext/assets3d/al_Asset.hpp"
#include <algorithm>
#include <cstdint>
#include <vector>
#include "soundSynth.cpp"
#include "al/app/al_GUIDomain.hpp"
#include "al/sphere/al_SphereUtils.hpp"

using namespace gam;
using namespace al;
using namespace std;
#define FFT_SIZE 4048
static bool fullscreen = true;

class BrainApp : public App
{
public:
  SynthGUIManager<FMWT> synthManager{"synth"};
  int midiNote;
  FMWT fm;
  Mesh mSpectrogram;
  vector<float> spectrum;
  // bool showGUI = true;
  // bool half_draw = false;
  // bool diff_color = false;
  // bool navi = false;
  bool showSpectro = false;
  gam::STFT stft = gam::STFT(FFT_SIZE, FFT_SIZE / 4, 0, gam::HANN, gam::MAG_FREQ);
  Scene *ascene1{nullptr};
  Scene *ascene2{nullptr};
  Vec3f scene_min, scene_max, scene_center;
  Vec3f scene_min_half, scene_max_half, scene_center_half;
  vector<Mesh> mMesh; 
  vector<Mesh> mMesh_half;
  Mesh center_of_part_mesh;

  Texture transparent_tex;
  double time;
  Light light;
  // ControlGUI *gui;

  ParameterBool showGUI{"showGUI", "", 1.0};
  ParameterBool half_draw{"Show Half Brain", "", 0.0};
  ParameterBool diff_color{"Differentiate Parts (color)", "", 0.0};
  ParameterBool navi{"Navigate", "", 1.0};
  ParameterInt upper_draw{"Draw Control", "",204,"", 1, 204};

  void PlatformSetupSize()
  {
    int total_width, total_height;
    al::sphere::getFullscreenDimension(&total_width, &total_height);
    std::cout << total_width << ", " << total_height << std::endl;
    dimensions(0, 0, total_width, total_height);
  }
  virtual void onInit() override
  {
    if (al_get_hostname() == "moxi" || fullscreen)
    {
      PlatformSetupSize();
    }    
    nav().pos(0, 3, 30);
    lens().near(0.1).far(1000).fovy(60);

    // Generate a texture with an alpha channel for transparency
    transparent_tex.create2D(256, 256, Texture::RGBA8);
    int Nx = transparent_tex.width();
    int Ny = transparent_tex.height();
    vector<Colori> pix;
    pix.resize(Nx * Ny);

    // Import .obj file for the mesh
    // std::string fileName = "../obj/ducky.obj";
    std::string fileName = "../obj/brain.obj";
    ascene1 = Scene::import(fileName);
    if (!ascene1)
    {
      printf("error reading %s\n", fileName.c_str());
      return;
    }
    ascene1->getBounds(scene_min, scene_max);
    scene_center = (scene_min + scene_max) / 2.f;
    // ascene->print();
    mMesh.resize(ascene1->meshes());
    for (int i = 0; i < ascene1->meshes(); i += 1)
    {
      ascene1->mesh(i, mMesh[i]);
    }
    transparent_tex.submit(mMesh.data());

    // Half
    fileName = "../obj/brain_half.obj";
    ascene2 = Scene::import(fileName);
    if (!ascene2)
    {
      printf("error reading %s\n", fileName.c_str());
      return;
    }
    ascene2->getBounds(scene_min_half, scene_max_half);
    scene_center_half = (scene_min_half + scene_max_half) / 2.f;
    // ascene->print();
    mMesh_half.resize(ascene2->meshes());
    for (int i = 0; i < ascene2->meshes(); i += 1)
    {
      ascene2->mesh(i, mMesh_half[i]);
    }
    transparent_tex.submit(mMesh.data());

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
    // Initialize GUI and Parameter callbacks
    // if (isPrimary()) 
    {
      // auto guiDomain = GUIDomain::enableGUI(defaultWindowDomain());
      // gui = &guiDomain->newGUI();
      // *gui << showGUI << half_draw << diff_color << navi;
      // *gui << sequencer << recorder;
    }
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
    time += dt;
    Vec3f point_you_want_to_see = Vec3f(0,0,0);
    light.pos(nav().pos().x, nav().pos().y, nav().pos().z);
    light.dir(point_you_want_to_see.x, point_you_want_to_see.y, point_you_want_to_see.z);
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
    {
      g.pushMatrix();
      float tmp = std::max(scene_max[1] - scene_min[1], tmp);
      tmp = std::max(scene_max[2] - scene_min[2], tmp);
      tmp = 2.f / tmp;
      g.scale(tmp);
      g.scale(10);
      g.meshColor();
      float color_index;
      // g.texture();
      // transparent_tex.bind();
      if (half_draw)
      {
        g.translate(0,0.011,0);
        for (auto &m : mMesh_half)
        // half 137 parts of brain
        {
          if(color_index < upper_draw-67)
          {   
            color_index +=1;
            g.translate( 0.000001 * sin(color_index* time*0.1), 0.000001 * sin(time*0.2), 0.000001 * cos(color_index* time*0.8) );
            // g.color(RGB( (color_index+200) / 255., 174 / 255., 177 / 255.)); 
            if(diff_color){
              g.color(HSV( (color_index) / 204., (sin(time*0.05*color_index)*0.1*color_index + 157) / 255., 177 / 255.)); 
            } else{
              g.color(RGB( (0.1*color_index+230) / 255., 174 / 255., 177 / 255.)); 
            }
            // g.color(RGB(242 / 255., 174 / 255., 177 / 255.));
            g.draw(m);
          }
        }
      }
      else
      {
        g.pushMatrix();      
        for (auto &m : mMesh)
        // 204 parts of brain
        {
          if(color_index < upper_draw)
          {       
            color_index +=1;
          g.translate( 0.000001 * sin(color_index* time*0.1), 0.000001 * sin(time*0.2), 0.000001 * cos(color_index* time*0.8) );
            // g.color(RGB( (color_index+200) / 255., 174 / 255., 177 / 255.)); 
            if(diff_color){
              g.color(HSV( (color_index) / 204., (sin(time*0.05*color_index)*0.1*color_index + 157) / 255., 177 / 255.)); 
            } else{
              g.color(RGB( (0.1*color_index+230) / 255., 177/ 255., 177 / 255.)); 
            }
            // g.color(RGB(242 / 255., 174 / 255., 177 / 255.));
            g.draw(m);
            g.draw(center_of_part_mesh);
          }
        }
        g.popMatrix();

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
    case '0':
      diff_color = !diff_color;
      break;
    case '7':
      upper_draw = 204;
    case '8':
      upper_draw = upper_draw-1;
      if(upper_draw<2)
        upper_draw = 1;
      break;
    case '9':
      upper_draw = upper_draw+1;
      if(upper_draw>204)
        upper_draw = 204;
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
  BrainApp app;

  // Set up audiom
  app.configureAudio(48000., 512, 2, 0);

  app.start();
}
