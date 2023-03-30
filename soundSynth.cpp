using namespace gam;
using namespace al;
using namespace std;
#define FFT_SIZE 4048
#define BLOCK_SIZE 512
#define CHANNEL_COUNT 2
// tables for oscillator
gam::ArrayPow2<float> tbSaw(2048), tbSqr(2048), tbImp(2048), tbSin(2048),
    tbPls(2048), tb__1(2048), tb__2(2048), tb__3(2048), tb__4(2048);

class FMWT : public SynthVoice
{
public:
  // Unit generators
  gam::Pan<> mPan;
  gam::ADSR<> mAmpEnv;
  gam::ADSR<> mModEnv;
  gam::EnvFollow<> mEnvFollow;
  gam::ADSR<> mVibEnv;
  Reverb<float> reverb;
  gam::Biquad<> mFilter{};
  gam::STFT stft = gam::STFT(FFT_SIZE, FFT_SIZE / 4, 0, gam::HANN, gam::MAG_FREQ);
  // This time, let's use spectrograms for each notes as the visual components.
  Mesh mSpectrogram;
  vector<float> spectrum;
  Mesh waveformMesh;
  gam::Sine<> mod, mVib; // carrier, modulator sine oscillators
  gam::Osc<> car;
  double a = 0;
  double b = 0;
  double timepose = 10;
  Vec3f init_pose, init_rot;
  float init_rot_quant;
  // Additional members
  float mVibFrq;
  float mVibDepth;
  float mVibRise;
  int mtable;
  static const int numb_waveform = 9;
  Mesh mMesh[numb_waveform];
  bool wireframe = false;
  bool vertexLight = false;
  float waveformData[BLOCK_SIZE * CHANNEL_COUNT]{0};
  void init() override
  {
    reverb.bandwidth(0.6f); // Low-pass amount on input, in [0,1]
    reverb.damping(0.5f);   // High-frequency damping, in [0,1]
    reverb.decay(0.6f);     // Tail decay factor, in [0,1]
    // mSpectrogram.primitive(Mesh::LINE_STRIP);
    mSpectrogram.primitive(Mesh::POINTS);
    // waveformMesh.primitive(Mesh::LINE_STRIP);
    waveformMesh.primitive(Mesh::POINTS);
    // Diffusion amounts
    // Values near 0.7 are recommended. Moving further away from 0.7 will lead
    // to more distinct echoes.
    reverb.diffusion(0.76, 0.666, 0.707, 0.571);
    mModEnv.levels(0, 1, 1, 0);
    mVibEnv.levels(0, 1, 1, 0);
    mAmpEnv.sustainPoint(2);
    mFilter.type(gam::FilterType(0));
    spectrum.resize(FFT_SIZE / 2 + 1);

    createInternalTriggerParameter("frequency", 440, 10, 4000.0);
    createInternalTriggerParameter("amplitude", 0.1, 0.0, 1.0);
    createInternalTriggerParameter("attackTime", 0.1, 0.01, 3.0);
    createInternalTriggerParameter("releaseTime", 5., 0.1, 10.0);
    createInternalTriggerParameter("sustain", 0.65, 0.1, 1.0);
    // FM index
    createInternalTriggerParameter("idx1", 0.01, 0.0, 10.0);
    createInternalTriggerParameter("idx2", 7, 0.0, 10.0);
    createInternalTriggerParameter("idx3", 5, 0.0, 10.0);

    createInternalTriggerParameter("carMul", 1.7, 0.0, 20.0);
    createInternalTriggerParameter("modMul", 1.30, 0.0, 20.0);

    createInternalTriggerParameter("vibRate1", 0.01, 0.0, 10.0);
    createInternalTriggerParameter("vibRate2", 0.5, 0.0, 10.0);
    createInternalTriggerParameter("vibRise", 0, 0.0, 10.0);
    createInternalTriggerParameter("vibDepth", 0, 0.0, 10.0);

    createInternalTriggerParameter("pan", 0.0, -1.0, 1.0);
    createInternalTriggerParameter("table", 1, 0, 8);
    createInternalTriggerParameter("reverb", 0.95, 0., 1);
    // 	LOW_PASS,	HIGH_PASS
    createInternalTriggerParameter("filterType", 0, 0, 2); 
    createInternalTriggerParameter("filterFrequency", 1200, 5, 5000);
    createInternalTriggerParameter("filterResonance", 1.0, 0.1f, 20);

    // Table & Visual meshes
    // Now We have the mesh according to the waveform
    gam::addSinesPow<1>(tbSaw, 9, 1);
    addCone(mMesh[0],1, Vec3f(0,0,5), 40, 1); //tbSaw

    gam::addSinesPow<1>(tbSqr, 9, 2);
    addCube(mMesh[1]);  // tbSquare

    gam::addSinesPow<0>(tbImp, 9, 1);
    addPrism(mMesh[2],1,1,1,100); // tbImp

    gam::addSine(tbSin);
    addSphere(mMesh[3], 0.3, 16, 100); // tbSin

// About: addSines (dst, amps, cycs, numh)
// \param[out] dst		destination array
// \param[in] amps		harmonic amplitudes of series, size must be numh - A[]
// \param[in] cycs		harmonic numbers of series, size must be numh - C[]
// \param[in] numh		total number of harmonics
    float scaler = 0.15;
    float hscaler = 1;

    { //tbPls
      float A[] = {1, 1, 1, 1, 0.7, 0.5, 0.3, 0.1};
      gam::addSines(tbPls, A, 8); 
      addWireBox(mMesh[4],2);    // tbPls
    }
    { // tb__1 
      float A[] = {1, 0.4, 0.65, 0.3, 0.18, 0.08, 0, 0};
      float C[] = {1, 4, 7, 11, 15, 18, 0, 0 };
      gam::addSines(tb__1, A, C, 6);
      for (int i = 0; i < 7; i++){
        addWireBox(mMesh[5], scaler * A[i]*C[i], scaler * A[i+1]*C[i+1], 1 + 0.3*i);

        // addSphere(mMesh[5],scaler * A[i], 16, 30); // tb__1
      }
    }
    { // inharmonic partials
      float A[] = {0.5, 0.8, 0.7, 1, 0.3, 0.4, 0.2, 0.12};
      float C[] = {3, 4, 7, 8, 11, 12, 15, 16}; 
      gam::addSines(tb__2, A, C, 8); // tb__2
      for (int i = 0; i < 7; i++){
        addWireBox(mMesh[6], scaler * A[i]*C[i], scaler * A[i+1]*C[i+1], 1 + 0.3*i);
      }
    }
    { // inharmonic partials
      float A[] = {1, 0.7, 0.45, 0.3, 0.15, 0.08, 0 , 0};
      float C[] = {10, 27, 54, 81, 108, 135, 0, 0};
      gam::addSines(tb__3, A, C, 6); // tb__3
      for (int i = 0; i < 7; i++){
        addWireBox(mMesh[7], scaler * A[i]*C[i], scaler * A[i+1]*C[i+1], 1 + 0.3*i);
      }
    }
  { // harmonics 20-27
      float A[] = {0.2, 0.4, 0.6, 1, 0.7, 0.5, 0.3, 0.1};
      gam::addSines(tb__4, A, 8, 20); // tb__4
      for (int i = 0; i < 7; i++){
        addWireBox(mMesh[8], hscaler * A[i], hscaler * A[i+1], 1 + 0.3*i);
      }
    }

    // Scale and generate normals
    for (int i = 0; i < numb_waveform; ++i) {
      mMesh[i].scale(0.4);

      int Nv = mMesh[i].vertices().size();
      for (int k = 0; k < Nv; ++k) {
        mMesh[i].color(HSV(float(k) / Nv, 0.3, 1));
      }

      if (!vertexLight && mMesh[i].primitive() == Mesh::TRIANGLES) {
        mMesh[i].decompress();
      }
      mMesh[i].generateNormals();
    }
  }

  //
  void onProcess(AudioIOData &io) override
  {
    mVib.freq(mVibEnv());
    float carBaseFreq =
        getInternalParameterValue("frequency") * getInternalParameterValue("carMul");
    float modScale = getInternalParameterValue("frequency") * getInternalParameterValue("modMul");
    float amp = getInternalParameterValue("amplitude") * 0.01;
    while (io())
    {
      mVib.freq(mVibEnv());
      car.freq((1 + mVib() * mVibDepth) * carBaseFreq +
               mod() * mModEnv() * modScale);
      float s1 = car() * mAmpEnv() * amp;
      // Compute two wet channels of reverberation
      float wet1, wet2;
      reverb(s1, wet1, wet2);
      float filtered, out1, out2;
      filtered = mFilter(wet1);

      mEnvFollow(filtered);
      mPan(filtered, out1, out2);
      io.out(0) += out1;
      io.out(1) += out2;
      // STFT for each notes
      if (stft(s1))
      { // Loop through all the frequency bins
          for (unsigned k = 0; k < stft.numBins(); ++k)
          {
              // Here we simply scale the complex sample
              spectrum[k] = tanh(pow(stft.bin(k).real(), 1.3));
          }
      }
      // waveform 
      memcpy(&waveformData, io.outBuffer(), BLOCK_SIZE * CHANNEL_COUNT * sizeof(float));
    }
    if (mAmpEnv.done() && (mEnvFollow.value() < 0.001))
      free();
  }

  void onProcess(Graphics &g) override
  {
    a += 0.29;
    b += 0.23;
    timepose -= 0.6;
    int shape = getInternalParameterValue("table");
    g.polygonMode(wireframe ? GL_LINE : GL_FILL);
    // light.pos(0, 0, 0);
    gl::depthTesting(true);
    mSpectrogram.reset();
    // mSpectrogram.primitive(Mesh::LINE_STRIP);

    for (int i = 0; i < FFT_SIZE / 2; i++)
    {
        mSpectrogram.color(HSV(spectrum[i] * 1000 + al::rnd::uniform()));
        mSpectrogram.vertex(i, spectrum[i], 0.0);
    }
    // wavefrom
    float chRatio = 1 / float(CHANNEL_COUNT);
    float hSegment = chRatio * 10;

    for(int ch = 0; ch < CHANNEL_COUNT; ch++) {
      waveformMesh.reset();
      float yBase = (CHANNEL_COUNT - ch - 1) * hSegment;
      for(int i = 0; i < BLOCK_SIZE; i++) {
        float x = 10 * (i / float(BLOCK_SIZE));
        float y = hSegment * ((waveformData[ch * BLOCK_SIZE + i] + 1)/ 2.0) + yBase;
        
        waveformMesh.vertex(x, y);
      }
    }
    g.pushMatrix();
    g.depthTesting(true);
    g.lighting(true);
    g.translate(init_pose);
    // g.translate(getInternalParameterValue("freq") / 100, getInternalParameterValue("freq") / 100 + 2, 0);
    g.rotate(init_rot_quant, init_rot);
    g.rotate(mVib(), Vec3f(0, 1, 0));
    g.rotate(mVib() * mVibDepth, Vec3f(1));
    float scaling = getInternalParameterValue("amplitude");
    // g.scale(10, 0.1,0.1);
    // g.color(HSV( (getInternalParameterValue("modMul") + init_rot_quant) / 360, 100 * mEnvFollow.value(), mEnvFollow.value()+0.9));
    // g.draw(mMesh[shape]);
    g.pointSize(3);
    // g.scale(10.0 / FFT_SIZE, 500, 1.0);
    g.color(HSV( (getInternalParameterValue("modMul") + init_rot_quant) / 360, 0.8 +100 * mEnvFollow.value(), 100.));
    // g.draw(mSpectrogram);
    g.scale(0.3, 1, 1);

    g.draw(waveformMesh);
    g.popMatrix();
  }

  void onTriggerOn() override
  {
    timepose = 10;
    init_pose = Vec3f(al::rnd::uniform(-3., 3.), al::rnd::uniform(-2., 4.), al::rnd::uniform(-5., 5.));
    init_rot = Vec3f(al::rnd::uniform(-5., 5.), al::rnd::uniform(-2., 5.), al::rnd::uniform(-5., 5.));
    init_rot_quant = al::rnd::uniform(180., 270.);
    mAmpEnv.reset();
    mVibEnv.reset();
    mModEnv.reset();
    mVib.phase(0);
    mod.phase(0);
    reverb.zero();
    mFilter.zero();
    updateFromParameters();
    updateWaveform();

    float modFreq =
        getInternalParameterValue("freq") * getInternalParameterValue("modMul");
    mod.freq(modFreq);
  }
  void onTriggerOff() override
  {
    mAmpEnv.triggerRelease();
    mModEnv.triggerRelease();
    mVibEnv.triggerRelease();
  }

  void updateFromParameters()
  {
    mFilter.type(gam::FilterType(getInternalParameterValue("filterType")));
    mFilter.freq(getInternalParameterValue("filterFrequency"));
    mFilter.res(getInternalParameterValue("filterResonance"));
    mModEnv.levels()[0] = getInternalParameterValue("idx1");
    mModEnv.levels()[1] = getInternalParameterValue("idx2");
    mModEnv.levels()[2] = getInternalParameterValue("idx2");
    mModEnv.levels()[3] = getInternalParameterValue("idx3");

    mAmpEnv.attack(getInternalParameterValue("attackTime"));
    mAmpEnv.release(getInternalParameterValue("releaseTime"));
    mAmpEnv.sustain(getInternalParameterValue("sustain"));

    mModEnv.lengths()[0] = getInternalParameterValue("attackTime");
    mModEnv.lengths()[3] = getInternalParameterValue("releaseTime");

    mVibEnv.levels(getInternalParameterValue("vibRate1"),
                   getInternalParameterValue("vibRate2"),
                   getInternalParameterValue("vibRate2"),
                   getInternalParameterValue("vibRate1"));
    mVibEnv.lengths()[0] = getInternalParameterValue("vibRise");
    mVibEnv.lengths()[1] = getInternalParameterValue("vibRise");
    mVibEnv.lengths()[3] = getInternalParameterValue("vibRise");
    mVibDepth = getInternalParameterValue("vibDepth");
    
    mPan.pos(getInternalParameterValue("pan"));
    reverb.decay(getInternalParameterValue("reverb"));
  }
  void updateWaveform(){
        // Map table number to table in memory
    switch (int(getInternalParameterValue("table"))) {
      case 0:
        car.source(tbSaw);
        break;
      case 1:
        car.source(tbSqr);
        break;
      case 2:
        car.source(tbImp);
        break;
      case 3:
        car.source(tbSin);
        break;
      case 4:
        car.source(tbPls);
        break;
      case 5:
        car.source(tb__1);
        break;
      case 6:
        car.source(tb__2);
        break;
      case 7:
        car.source(tb__3);
        break;
      case 8:
        car.source(tb__4);
        break;
    }
  }
};
