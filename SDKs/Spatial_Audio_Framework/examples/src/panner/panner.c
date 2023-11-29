/*
 * Copyright 2017-2018 Leo McCormack
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
 * INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 * LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
 * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 * PERFORMANCE OF THIS SOFTWARE.
 */

/**
 * @file panner.c
 * @brief A frequency-dependent 3D panner based on the Vector-base Amplitude
 *        Panning (VBAP) method [1], with an optional spread control [2].
 *
 * Depending on the listening room, it may be beneficial to employ amplitude-
 * normalised gains for low frequencies, and energy-normalised gains for high
 * frequencies. Therefore, this VBAP implementation also uses the method
 * described in [3], to do just that.
 *
 * @see [1] Pulkki, V. (1997). Virtual sound source positioning using vector
 *          base amplitude panning. Journal of the audio engineering society,
 *          45(6), 456-466.
 * @see [2] Pulkki, V. (1999). Uniform spreading of amplitude panned virtual
 *          sources. In Proceedings of the 1999 IEEE Workshop on Applications of
 *          Signal Processing to Audio and Acoustics. WASPAA'99 (Cat. No.
 *          99TH8452) (pp. 187-190). IEEE.
 * @see [3] Laitinen, M., Vilkamo, J., Jussila, K., Politis, A., Pulkki, V.
 *          (2014). Gain normalisation in amplitude panning as a function of
 *          frequency and room reverberance. 55th International Conference of
 *          the AES. Helsinki, Finland.
 *
 * @author Leo McCormack
 * @date 25.09.2017
 * @license ISC
 */
 
#include "panner_internal.h"
void panner_create
(
    void ** const phPan
)
{
    panner_data* pData = (panner_data*)malloc1d(sizeof(panner_data));
    *phPan = (void*)pData;
    int ch, dummy;

    /* default user parameters */
    panner_loadSourcePreset(SOURCE_CONFIG_PRESET_DEFAULT, pData->src_dirs_deg, &(pData->new_nSources), &(dummy)); /*check setStateInformation if you change default preset*/
    pData->nSources = pData->new_nSources;
    
    
    /* time-frequency transform + buffers */
    pData->hSTFT = NULL;

    int n, i, band;

    /* Default user parameters */
    pData->masterOrder = pData->new_masterOrder = SH_ORDER_FIRST;
    for (band = 0; band < HYBRID_BANDS; band++) {
        pData->analysisOrderPerBand[band] = pData->masterOrder;
        pData->pmapEQ[band] = 1.0f;
    }
    pData->covAvgCoeff = 0.3f;
    pData->pmapAvgCoeff = 0.0f;
    pData->nSources = 1;
    pData->HFOVoption = HFOV_360;
    pData->aspectRatioOption = ASPECT_RATIO_2_1;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D;

    pData->spPWD = NULL;

    afSTFT_create(&(pData->hSTFT), MAX_NUM_SH_SIGNALS, 0, HOP_SIZE, 0, 1, AFSTFT_BANDS_CH_TIME);
    pData->SHframeTD = (float**)malloc2d(MAX_NUM_SH_SIGNALS, PANNER_FRAME_SIZE, sizeof(float));
    pData->SHframeTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_SH_SIGNALS, TIME_SLOTS, sizeof(float_complex));


    /* codec data */
    pData->pars = (powermap_codecPars*)malloc1d(sizeof(powermap_codecPars));
    powermap_codecPars* pars = pData->pars;
    pars->interp_dirs_deg = NULL;
    for (n = 0; n < MAX_SH_ORDER; n++) {
        pars->Y_grid[n] = NULL;
        pars->Y_grid_cmplx[n] = NULL;
    }
    pars->interp_table = NULL;

    /* internal */
    pData->progressBar0_1 = 0.0f;
    pData->progressBarText = malloc1d(PROGRESSBARTEXT_CHAR_LENGTH * sizeof(char));
    strcpy(pData->progressBarText, "");
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
    pData->procStatus = PROC_STATUS_NOT_ONGOING;
    pData->dispWidth = 140;

    /* display */
    pData->pmap = NULL;
    pData->prev_pmap = NULL;
    for (i = 0; i < NUM_DISP_SLOTS; i++)
        pData->pmap_grid[i] = NULL;
    pData->pmapReady = 0;
    pData->recalcPmap = 1;


}

void panner_destroy
(
    void ** const phPan
)
{
    panner_data* pData = (panner_data*)(*phPan);
    powermap_codecPars* pars;
    int i;

    if (pData != NULL) {
        /* not safe to free memory during intialisation/processing loop */
        while (pData->codecStatus == CODEC_STATUS_INITIALISING ||
            pData->procStatus == PROC_STATUS_ONGOING) {
            SAF_SLEEP(10);
        }
        sphPWD_destroy(&(pData->spPWD));

        pars = pData->pars;

        /* free afSTFT and buffers */
        afSTFT_destroy(&(pData->hSTFT));
        free(pData->SHframeTD);
        free(pData->SHframeTF);

        free(pData->pmap);
        free(pData->prev_pmap);
        for (i = 0; i < NUM_DISP_SLOTS; i++)
            free(pData->pmap_grid[i]);
        free(pars->interp_dirs_deg);
        for (i = 0; i < MAX_SH_ORDER; i++) {
            free(pars->Y_grid[i]);
            free(pars->Y_grid_cmplx[i]);
        }
        //sphPWD_destroy(&(pData->spPWD));
        free(pars->interp_table);
        free(pData->pars);
        free(pData->progressBarText);
        free(pData);
        pData = NULL;
    }
}

void panner_init
(
    void * const hPan,
    int          sampleRate
)
{
    panner_data *pData = (panner_data*)(hPan);
    powermap_codecPars* pars = pData->pars;
    
    /* define frequency vector */
    pData->fs = sampleRate;
    afSTFT_getCentreFreqs(pData->hSTFT, (float)sampleRate, HYBRID_BANDS, pData->freqVector);
    /* intialise parameters */
    memset(pData->Cx, 0, MAX_NUM_SH_SIGNALS * MAX_NUM_SH_SIGNALS * HYBRID_BANDS * sizeof(float_complex));
    if (pData->prev_pmap != NULL)
        memset(pData->prev_pmap, 0, pars->grid_nDirs * sizeof(float));
    pData->pmapReady = 0;
    pData->dispSlotIdx = 0;
}

void panner_initCodec
(
    void* const hPan
)
{
    panner_data *pData = (panner_data*)(hPan);
    
    if (pData->codecStatus != CODEC_STATUS_NOT_INITIALISED)
        return; /* re-init not required, or already happening */
    while (pData->procStatus == PROC_STATUS_ONGOING){
        /* re-init required, but we need to wait for the current processing loop to end */
        pData->codecStatus = CODEC_STATUS_INITIALISING; /* indicate that we want to init */
        SAF_SLEEP(10);
    }
    
    /* for progress bar */
    pData->codecStatus = CODEC_STATUS_INITIALISING;
    strcpy(pData->progressBarText,"Initialising");
    pData->progressBar0_1 = 0.0f;
    
    /* reinit TFT if needed */
    panner_initTFT(hPan);
    panner_initAna(hPan);

    
    /* done! */
    strcpy(pData->progressBarText,"Done!");
    pData->progressBar0_1 = 1.0f;
    pData->codecStatus = CODEC_STATUS_INITIALISED;
    
}

float panner_beamform(const float X[], const float Y[], const float Z[], void* const bPan, int index, int num_samples) {
    Direction* bData = (Direction*)(bPan);

    float energy = 0.0f;
    float ux = bData[index].cosPhi * bData[index].cosTheta;
    float uy = bData[index].cosPhi * bData[index].sinTheta;
    float uz = bData[index].sinPhi;
    for (int i = 0; i < num_samples; i++) {
        float sample_energy = ux * X[i] + uy * Y[i] + uz * Z[i];
        energy += sample_energy * sample_energy; // Assuming energy is sum of squared sample energies
    }
    if (energy == INFINITE || energy == INFINITY) {
        printf("chuj");
    }
    return energy;
}

void panner_beamformer_process(const float X[], const float Y[], const float Z[], int numSamples, int loudNum, void * const bPan, void * const hPan)
{
    Direction* bData = (Direction*)(bPan);
    bData->status = PROC_STATUS_ONGOING;
    panner_data* pData = (panner_data*)(hPan);
    float maxEnergy = -1e9f;
    float bestTheta = 0.0f, bestPhi = 0.0f;

//#pragma omp parallel for reduction(max:maxEnergy)
    for (int i = 0; i < 360*180; i++) {
        float energy = panner_beamform(X, Y, Z, bPan, i, numSamples);
        if (energy > maxEnergy) {
//#pragma omp critical
            {
                if (energy > maxEnergy) {
                    maxEnergy = energy;
                    bestTheta = bData[i].theta;
                    bestPhi = bData[i].phi;
                }
            }
        }
    }
    bestTheta = (bestTheta * 180.0f / PI) < 0 ? (bestTheta * 180.0f / PI) + 360.0f : (bestTheta * 180.0f / PI);
    bestPhi = bestPhi * 180.0f / PI;
    panner_setLoudspeakerAzi_deg(hPan, loudNum, bestTheta);
    panner_setLoudspeakerElev_deg(hPan, loudNum, bestPhi);
    bData->status = PROC_STATUS_NOT_ONGOING;
}



float toRadians(float degrees)
{
    return degrees * (PI / 180.0);
}

void calculateCoordinates(float distance, float azimuth, float* x, float* y) //we get x and y
{
    float azimuthRadians = toRadians(azimuth);
    *x = distance * cos(azimuthRadians);
    *y = distance * sin(azimuthRadians);
}


void panner_process(void* const hPan,
    const float* const* inputs,
    int                  nInputs,
    int                  nSamples,
    int                  loudNum,
    int                  isPlaying
)
{
    panner_data* pData = (panner_data*)(hPan);
    powermap_codecPars* pars = pData->pars;
    int s, i, j, ch, band, nSH_order, order_band, nSH_maxOrder, maxOrder;
    float C_grp_trace, pmapEQ_band;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    float_complex new_Cx[MAX_NUM_SH_SIGNALS * MAX_NUM_SH_SIGNALS];
    float_complex C_grp[MAX_NUM_SH_SIGNALS * MAX_NUM_SH_SIGNALS];

    /* local parameters */
    int analysisOrderPerBand[HYBRID_BANDS];
    int nSources, masterOrder, nSH;
    float covAvgCoeff, pmapAvgCoeff;
    float pmapEQ[HYBRID_BANDS];
    NORM_TYPES norm;
    CH_ORDER chOrdering;
    memcpy(analysisOrderPerBand, pData->analysisOrderPerBand, HYBRID_BANDS * sizeof(int));
    memcpy(pmapEQ, pData->pmapEQ, HYBRID_BANDS * sizeof(float));
    norm = pData->norm;
    chOrdering = pData->chOrdering;
    nSources = pData->nSources;
    covAvgCoeff = SAF_MIN(pData->covAvgCoeff, MAX_COV_AVG_COEFF);
    pmapAvgCoeff = pData->pmapAvgCoeff;
    masterOrder = pData->masterOrder;
    nSH = (masterOrder + 1) * (masterOrder + 1);

    pData->procStatus = PROC_STATUS_ONGOING;

    /* Load time-domain data */
    for (ch = 0; ch < nSH; ch++)
        memcpy(pData->SHframeTD[ch], inputs[ch], nSamples * sizeof(float));

    /* account for input channel order */
    switch (chOrdering) {
    case CH_ACN:  /* already ACN */ break; /* Otherwise, convert to ACN... */
    case CH_FUMA: convertHOAChannelConvention(FLATTEN2D(pData->SHframeTD), masterOrder, nSamples, HOA_CH_ORDER_FUMA, HOA_CH_ORDER_ACN); break;
    }

    /* account for input normalisation scheme */
    switch (norm) {
    case NORM_N3D:  /* already in N3D, do nothing */ break; /* Otherwise, convert to N3D... */
    case NORM_SN3D: convertHOANormConvention(FLATTEN2D(pData->SHframeTD), masterOrder, nSamples, HOA_NORM_SN3D, HOA_NORM_N3D); break;
    case NORM_FUMA: convertHOANormConvention(FLATTEN2D(pData->SHframeTD), masterOrder, nSamples, HOA_NORM_FUMA, HOA_NORM_N3D); break;
    }

    /* apply the time-frequency transform */
    afSTFT_forward_knownDimensions(pData->hSTFT, pData->SHframeTD, nSamples, MAX_NUM_SH_SIGNALS, TIME_SLOTS, pData->SHframeTF);

    /* Update covarience matrix per band */
    for (band = 0; band < HYBRID_BANDS; band++) {
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nSH, nSH, TIME_SLOTS, &calpha,
            FLATTEN2D(pData->SHframeTF[band]), TIME_SLOTS,
            FLATTEN2D(pData->SHframeTF[band]), TIME_SLOTS, &cbeta,
            new_Cx, nSH);

        /* average over time */
        cblas_sscal(nSH * nSH * 2, covAvgCoeff, (float*)pData->Cx[band], 1);
        cblas_saxpy(nSH * nSH * 2, 1.0f - covAvgCoeff, (float*)new_Cx, 1, (float*)pData->Cx[band], 1);
    }

    /* determine maximum analysis order */
    maxOrder = 1;
    for (i = 0; i < HYBRID_BANDS; i++)
        maxOrder = SAF_MAX(maxOrder, SAF_MIN(analysisOrderPerBand[i], masterOrder));
    nSH_maxOrder = (maxOrder + 1) * (maxOrder + 1);

    /* group covarience matrices */
    memset(C_grp, 0, nSH_maxOrder * nSH_maxOrder * sizeof(float_complex));
    for (band = 0; band < HYBRID_BANDS; band++) {
        order_band = SAF_MAX(SAF_MIN(pData->analysisOrderPerBand[band], masterOrder), 1);
        nSH_order = (order_band + 1) * (order_band + 1);
        pmapEQ_band = SAF_MIN(SAF_MAX(pmapEQ[band], 0.0f), 2.0f);
        for (i = 0; i < nSH_order; i++)
            for (j = 0; j < nSH_order; j++)
                C_grp[i * nSH_maxOrder + j] = ccaddf(C_grp[i * nSH_maxOrder + j], crmulf(pData->Cx[band][i * nSH + j], 1e3f * pmapEQ_band));
    }

    /* generate powermap */
    C_grp_trace = 0.0f;
    for (i = 0; i < nSH_maxOrder; i++)
        C_grp_trace += crealf(C_grp[i * nSH_maxOrder + i]);
    int indd;
    sphPWD_compute(pData->spPWD, (float_complex*)C_grp, 1, NULL, &indd);

    panner_setLoudspeakerAzi_deg(hPan, loudNum, pars->interp_dirs_deg[indd*2] > 0 ? pars->interp_dirs_deg[indd * 2] : pars->interp_dirs_deg[indd * 2] + 360);
    panner_setLoudspeakerElev_deg(hPan, loudNum, pars->interp_dirs_deg[indd*2 + 1]);
      
    pData->procStatus = PROC_STATUS_NOT_ONGOING;
}

void panner_requestPmapUpdate(void* const hPan)
{
    panner_data* pData = (panner_data*)(hPan);
    pData->recalcPmap = 1;
}

int panner_get_directions(void* const hPan) {
    panner_data* pData = (panner_data*)(hPan);
    return pData->pmapReady;
}
//void panner_process
//(
//    void        *  const hPan,
//    const float *const * inputs,
//    float       ** const outputs,
//    int                  nInputs,
//    int                  nOutputs,
//    int                  nSamples
//)
//{
//    panner_data *pData = (panner_data*)(hPan);
//    int t, ch, ls, i, band, nSources, nLoudspeakers, N_azi, aziIndex, elevIndex, idx3d, idx2D;
//    float aziRes, elevRes, pv_f, gains3D_sum_pvf, gains2D_sum_pvf, Rxyz[3][3], hypotxy;
//    float src_dirs[MAX_NUM_INPUTS][2], pValue[HYBRID_BANDS], gains3D[MAX_NUM_OUTPUTS], gains2D[MAX_NUM_OUTPUTS];
//    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
//    float_complex outputTemp[MAX_NUM_OUTPUTS][TIME_SLOTS];
//
//    /* copy user parameters to local variables */
//    memcpy(src_dirs, pData->src_dirs_deg, MAX_NUM_INPUTS*2*sizeof(float));
//    memcpy(pValue, pData->pValue, HYBRID_BANDS*sizeof(float));
//    nSources = pData->nSources;
//    nLoudspeakers = pData->nLoudpkrs;
//
//    /* apply panner */
//    if ((nSamples == PANNER_FRAME_SIZE) && (pData->vbap_gtable != NULL) && (pData->codecStatus == CODEC_STATUS_INITIALISED) ) {
//        pData->procStatus = PROC_STATUS_ONGOING;
//
//        /* Load time-domain data */
//        for(i=0; i < SAF_MIN(nSources,nInputs); i++)
//            utility_svvcopy(inputs[i], PANNER_FRAME_SIZE, pData->inputFrameTD[i]);
//        for(; i<MAX_NUM_INPUTS; i++)
//            memset(pData->inputFrameTD[i], 0, PANNER_FRAME_SIZE * sizeof(float));
//
//        /* Apply time-frequency transform (TFT) */
//        afSTFT_forward_knownDimensions(pData->hSTFT, pData->inputFrameTD, PANNER_FRAME_SIZE, MAX_NUM_INPUTS, TIME_SLOTS, pData->inputframeTF);
//        memset(FLATTEN3D(pData->outputframeTF), 0, HYBRID_BANDS*MAX_NUM_OUTPUTS*TIME_SLOTS * sizeof(float_complex));
//        memset(outputTemp, 0, MAX_NUM_OUTPUTS*TIME_SLOTS * sizeof(float_complex));
//
//        /* Rotate source directions */
//        if(pData->recalc_M_rotFLAG){
//            yawPitchRoll2Rzyx (pData->yaw, pData->pitch, pData->roll, 0, Rxyz);
//            for(i=0; i<nSources; i++){
//                pData->src_dirs_xyz[i][0] = cosf(DEG2RAD(pData->src_dirs_deg[i][1])) * cosf(DEG2RAD(pData->src_dirs_deg[i][0]));
//                pData->src_dirs_xyz[i][1] = cosf(DEG2RAD(pData->src_dirs_deg[i][1])) * sinf(DEG2RAD(pData->src_dirs_deg[i][0]));
//                pData->src_dirs_xyz[i][2] = sinf(DEG2RAD(pData->src_dirs_deg[i][1]));
//                pData->recalc_gainsFLAG[i] = 1;
//            }
//            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSources, 3, 3, 1.0f,
//                        (float*)(pData->src_dirs_xyz), 3,
//                        (float*)Rxyz, 3, 0.0f,
//                        (float*)(pData->src_dirs_rot_xyz), 3);
//            for(i=0; i<nSources; i++){
//                hypotxy = sqrtf(powf(pData->src_dirs_rot_xyz[i][0], 2.0f) + powf(pData->src_dirs_rot_xyz[i][1], 2.0f));
//                pData->src_dirs_rot_deg[i][0] = RAD2DEG(atan2f(pData->src_dirs_rot_xyz[i][1], pData->src_dirs_rot_xyz[i][0]));
//                pData->src_dirs_rot_deg[i][1] = RAD2DEG(atan2f(pData->src_dirs_rot_xyz[i][2], hypotxy));
//            }
//            pData->recalc_M_rotFLAG = 0;
//        }
//
//        /* Apply VBAP Panning */
//        if(pData->output_nDims == 3){/* 3-D case */
//            aziRes = (float)pData->vbapTableRes[0];
//            elevRes = (float)pData->vbapTableRes[1];
//            N_azi = (int)(360.0f / aziRes + 0.5f) + 1;
//            for (ch = 0; ch < nSources; ch++) {
//                /* recalculate frequency dependent panning gains */
//                if(pData->recalc_gainsFLAG[ch]){
//                    //aziIndex = (int)(matlab_fmodf(pData->src_dirs_deg[ch][0] + 180.0f, 360.0f) / aziRes + 0.5f);
//                    //elevIndex = (int)((pData->src_dirs_deg[ch][1] + 90.0f) / elevRes + 0.5f);
//                    aziIndex = (int)(matlab_fmodf(pData->src_dirs_rot_deg[ch][0] + 180.0f, 360.0f) / aziRes + 0.5f);
//                    elevIndex = (int)((pData->src_dirs_rot_deg[ch][1] + 90.0f) / elevRes + 0.5f);
//                    idx3d = elevIndex * N_azi + aziIndex;
//                    for (ls = 0; ls < nLoudspeakers; ls++)
//                        gains3D[ls] =  pData->vbap_gtable[idx3d*nLoudspeakers+ls];
//                    for (band = 0; band < HYBRID_BANDS; band++){
//                        /* apply pValue per frequency */
//                        pv_f = pData->pValue[band];
//                        if(pv_f != 2.0f){
//                            gains3D_sum_pvf = 0.0f;
//                            for (ls = 0; ls < nLoudspeakers; ls++)
//                                gains3D_sum_pvf += powf(SAF_MAX(gains3D[ls], 0.0f), pv_f);
//                            gains3D_sum_pvf = powf(gains3D_sum_pvf, 1.0f/(pv_f+2.23e-9f));
//                            for (ls = 0; ls < nLoudspeakers; ls++)
//                                pData->G_src[band][ch][ls] = cmplxf(gains3D[ls] / (gains3D_sum_pvf+2.23e-9f), 0.0f);
//                        }
//                        else
//                            for (ls = 0; ls < nLoudspeakers; ls++)
//                                pData->G_src[band][ch][ls] = cmplxf(gains3D[ls], 0.0f);
//                    }
//                    pData->recalc_gainsFLAG[ch] = 0;
//                }
//            }
//            /* apply panning gains */
//            for (band = 0; band < HYBRID_BANDS; band++) {
//                cblas_cgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nLoudspeakers, TIME_SLOTS, nSources, &calpha,
//                    pData->G_src[band], MAX_NUM_OUTPUTS,
//                    FLATTEN2D(pData->inputframeTF[band]), TIME_SLOTS, &cbeta,
//                    outputTemp, TIME_SLOTS);
//                for (i = 0; i < nLoudspeakers; i++)
//                    for (t = 0; t < TIME_SLOTS; t++)
//                        pData->outputframeTF[band][i][t] = ccaddf(pData->outputframeTF[band][i][t], outputTemp[i][t]);
//            }
//        }
//        else{/* 2-D case */
//            aziRes = (float)pData->vbapTableRes[0];
//            for (ch = 0; ch < nSources; ch++) {
//                /* recalculate frequency dependent panning gains */
//                if(pData->recalc_gainsFLAG[ch]){
//                    //idx2D = (int)((matlab_fmodf(pData->src_dirs_deg[ch][0]+180.0f,360.0f)/aziRes)+0.5f);
//                    idx2D = (int)((matlab_fmodf(pData->src_dirs_rot_deg[ch][0]+180.0f,360.0f)/aziRes)+0.5f);
//                    for (ls = 0; ls < nLoudspeakers; ls++)
//                        gains2D[ls] = pData->vbap_gtable[idx2D*nLoudspeakers+ls];
//                    for (band = 0; band < HYBRID_BANDS; band++){
//                        /* apply pValue per frequency */
//                        pv_f = pData->pValue[band];
//                        if(pv_f != 2.0f){
//                            gains2D_sum_pvf = 0.0f;
//                            for (ls = 0; ls < nLoudspeakers; ls++)
//                                gains2D_sum_pvf += powf(SAF_MAX(gains2D[ls], 0.0f), pv_f);
//                            gains2D_sum_pvf = powf(gains2D_sum_pvf, 1.0f/(pv_f+2.23e-9f));
//                            for (ls = 0; ls < nLoudspeakers; ls++)
//                                pData->G_src[band][ch][ls] = cmplxf(gains2D[ls] / (gains2D_sum_pvf+2.23e-9f), 0.0f);
//                        }
//                        else
//                            for (ls = 0; ls < nLoudspeakers; ls++)
//                                pData->G_src[band][ch][ls] = cmplxf(gains2D[ls], 0.0f);
//                    }
//                    pData->recalc_gainsFLAG[ch] = 0;
//                }
//
//                /* apply panning gains */
//                for (band = 0; band < HYBRID_BANDS; band++){
//                    for (ls = 0; ls < nLoudspeakers; ls++)
//                        for (t = 0; t < TIME_SLOTS; t++)
//                            pData->outputframeTF[band][ls][t] = ccaddf(pData->outputframeTF[band][ls][t], ccmulf(pData->inputframeTF[band][ch][t], pData->G_src[band][ch][ls]));
//                }
//            }
//        }
//
//        /* scale by sqrt(number of sources) */
//        for (band = 0; band < HYBRID_BANDS; band++)
//            cblas_sscal(/*re+im*/2*nLoudspeakers*TIME_SLOTS, 1.0f/sqrtf((float)nSources), (float*)FLATTEN2D(pData->outputframeTF[band]), 1);
//
//        /* inverse-TFT and copy to output */
//        afSTFT_backward_knownDimensions(pData->hSTFT, pData->outputframeTF, PANNER_FRAME_SIZE, MAX_NUM_OUTPUTS, TIME_SLOTS, pData->outputFrameTD);
//        for (ch = 0; ch < SAF_MIN(nLoudspeakers, nOutputs); ch++)
//            utility_svvcopy(pData->outputFrameTD[ch], PANNER_FRAME_SIZE, outputs[ch]);
//        for (; ch < nOutputs; ch++)
//            memset(outputs[ch], 0, PANNER_FRAME_SIZE*sizeof(float));
//
//    }
//    else
//        for (ch=0; ch < nOutputs; ch++)
//            memset(outputs[ch],0, PANNER_FRAME_SIZE*sizeof(float));
//
//
//    pData->procStatus = PROC_STATUS_NOT_ONGOING;
//}


/* Set Functions */

void panner_refreshSettings(void* const hPan)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
    panner_setCodecStatus(hPan, CODEC_STATUS_NOT_INITIALISED);
}

void panner_setSourceAzi_deg(void* const hPan, int index, float newAzi_deg)
{
    panner_data *pData = (panner_data*)(hPan);
   // if(newAzi_deg>180.0f)
   //     newAzi_deg = -360.0f + newAzi_deg;
    newAzi_deg = SAF_MAX(newAzi_deg, 0.0f);
    newAzi_deg = SAF_MIN(newAzi_deg, 360.0f);
    if(pData->src_dirs_deg[index][0] != newAzi_deg){
        pData->src_dirs_deg[index][0] = newAzi_deg;
    }
}

void panner_setSourceElev_deg(void* const hPan, int index, float newElev_deg)
{
    panner_data *pData = (panner_data*)(hPan);
    newElev_deg = SAF_MAX(newElev_deg, -90.0f);
    newElev_deg = SAF_MIN(newElev_deg, 90.0f);
    if(pData->src_dirs_deg[index][1] != newElev_deg){
        pData->src_dirs_deg[index][1] = newElev_deg;
    }
}

void panner_setNumSources(void* const hPan, int new_nSources)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
    /* determine if TFT must be reinitialised */
    new_nSources = new_nSources > MAX_NUM_INPUTS ? MAX_NUM_INPUTS : new_nSources;
    if(pData->nSources != new_nSources){
        pData->new_nSources = new_nSources;
        panner_setCodecStatus(hPan, CODEC_STATUS_NOT_INITIALISED);
    }
}

void panner_setLoudspeakerAzi_deg(void* const hPan, int index, float newAzi_deg)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
   // if(newAzi_deg>180.0f)
  //      newAzi_deg = -360.0f + newAzi_deg;
    newAzi_deg = roundf(newAzi_deg * 10.0f) / 10.0f;
    newAzi_deg = SAF_MAX(newAzi_deg, 0.0f);
    newAzi_deg = SAF_MIN(newAzi_deg, 360.0f);
    if(pData->loudpkrs_dirs_deg[index][0] != newAzi_deg){
        pData->loudpkrs_dirs_deg[index][0] = newAzi_deg;
        panner_setCodecStatus(hPan, CODEC_STATUS_NOT_INITIALISED);
    }
}

void panner_setLoudspeakerElev_deg(void* const hPan, int index, float newElev_deg)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
    newElev_deg = roundf(newElev_deg * 10.0f) / 10.0f; //chanign to 0.0 XXXX
    newElev_deg = SAF_MAX(newElev_deg, -90.0f);
    newElev_deg = SAF_MIN(newElev_deg, 90.0f);
    if(pData->loudpkrs_dirs_deg[index][1] != newElev_deg){
        pData->loudpkrs_dirs_deg[index][1] = newElev_deg;
        panner_setCodecStatus(hPan, CODEC_STATUS_NOT_INITIALISED);
    }
}

void panner_setLoudspeakerDist_deg(void* const hPan, int index, float newDist) //XXXX
{
    panner_data* pData = (panner_data*)(hPan);
    int ch;
    newDist = roundf(newDist * 10.0f) / 10.0f;
    newDist = SAF_MAX(newDist, 0.0f); //min distance
    newDist = SAF_MIN(newDist, 10.0f); //max distance

    if (pData->loudpkrs_dirs_deg[index][2] != newDist) { //not sure if index is 2?
        pData->loudpkrs_dirs_deg[index][2] = newDist;
    }
}

void panner_setNumLoudspeakers(void* const hPan, int new_nLoudspeakers)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
    //for (int i = 0; i < 64; i++)
    //{
    //    for (int j = 0; i < 3; j++)
    //    {
    //        pData->loudpkrs_dirs_deg[i][0] = 0.0f;
    //    }
    //} //XXX don't work
    
    new_nLoudspeakers  = new_nLoudspeakers > MAX_NUM_OUTPUTS ? MAX_NUM_OUTPUTS : new_nLoudspeakers;
    if(pData->new_nLoudpkrs != new_nLoudspeakers){
        pData->new_nLoudpkrs = new_nLoudspeakers;
        panner_setCodecStatus(hPan, CODEC_STATUS_NOT_INITIALISED);
    }
}

void panner_setOutputConfigPreset(void* const hPan, int newPresetID)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch, dummy;
    panner_loadLoudspeakerPreset(newPresetID, pData->loudpkrs_dirs_deg, &(pData->new_nLoudpkrs), &dummy);
    panner_setCodecStatus(hPan, CODEC_STATUS_NOT_INITIALISED);
}

void panner_setInputConfigPreset(void* const hPan, int newPresetID)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch, dummy;
    panner_loadSourcePreset(newPresetID, pData->src_dirs_deg, &(pData->new_nSources), &dummy);
    panner_setCodecStatus(hPan, CODEC_STATUS_NOT_INITIALISED);
}


void panner_setChOrder(void* const hPan, int newOrder)
{
    panner_data* pData = (panner_data*)(hPan);
    if ((CH_ORDER)newOrder != CH_FUMA || pData->new_masterOrder == SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->chOrdering = (CH_ORDER)newOrder;
}

void panner_setNormType(void* const hPan, int newType)
{
    panner_data* pData = (panner_data*)(hPan);
    if ((NORM_TYPES)newType != NORM_FUMA || pData->new_masterOrder == SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->norm = (NORM_TYPES)newType;
}


/* Get Functions */

int panner_getFrameSize(void)
{
    return PANNER_FRAME_SIZE;
}

CODEC_STATUS panner_getCodecStatus(void* const hPan)
{
    panner_data *pData = (panner_data*)(hPan);
    return pData->codecStatus;
}

PROC_STATUS panner_getProcStatus(void* const hPan)
{
    panner_data* pData = (panner_data*)(hPan);
    return pData->procStatus;
}

PROC_STATUS panner_getBeamStatus(void* const bPan)
{
    Direction* bData = (Direction*)(bPan);
    return bData->status;
}

float panner_getProgressBar0_1(void* const hPan)
{
    panner_data *pData = (panner_data*)(hPan);
    return pData->progressBar0_1;
}

void panner_getProgressBarText(void* const hPan, char* text)
{
    panner_data *pData = (panner_data*)(hPan);
    memcpy(text, pData->progressBarText, PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
}

float panner_getSourceAzi_deg(void* const hPan, int index)
{
    panner_data *pData = (panner_data*)(hPan);
    return pData->src_dirs_deg[index][0];
}

float panner_getSourceElev_deg(void* const hPan, int index)
{
    panner_data *pData = (panner_data*)(hPan);
    return pData->src_dirs_deg[index][1];
}

int panner_getNumSources(void* const hPan)
{
    panner_data *pData = (panner_data*)(hPan);
    return pData->new_nSources;
}

int panner_getMaxNumSources()
{
    return MAX_NUM_INPUTS;
}

float panner_getLoudspeakerAzi_deg(void* const hPan, int index)
{
    panner_data *pData = (panner_data*)(hPan);
    pData->loudpkrs_dirs_deg[index][0] = roundf(pData->loudpkrs_dirs_deg[index][0] * 10.0f) / 10.0f;
    return pData->loudpkrs_dirs_deg[index][0];
}

float panner_getLoudspeakerElev_deg(void* const hPan, int index)
{
    panner_data *pData = (panner_data*)(hPan);
    pData->loudpkrs_dirs_deg[index][1] = roundf(pData->loudpkrs_dirs_deg[index][1] * 10.0f) / 10.0f;
    return pData->loudpkrs_dirs_deg[index][1];
}

float panner_getLoudspeakerDist_deg(void* const hPan, int index)
{
    panner_data* pData = (panner_data*)(hPan);
    pData->loudpkrs_dirs_deg[index][2] = roundf(pData->loudpkrs_dirs_deg[index][2] * 10.0f) / 10.0f;
	return pData->loudpkrs_dirs_deg[index][2];
}

int panner_getNumLoudspeakers(void* const hPan)
{
    panner_data* pData = (panner_data*)(hPan);

    return pData->new_nLoudpkrs;
}

int panner_getMaxNumLoudspeakers()
{
    return MAX_NUM_OUTPUTS;
}

int panner_getDAWsamplerate(void* const hPan)
{
    panner_data *pData = (panner_data*)(hPan);
    return pData->fs;
}


int panner_getChOrder(void* const hPan)
{
    panner_data* pData = (panner_data*)(hPan);
    return (int)pData->chOrdering;
}

int panner_getNormType(void* const hPan)
{
    panner_data* pData = (panner_data*)(hPan);
    return (int)pData->norm;
}

int panner_getProcessingDelay()
{
    return PANNER_FRAME_SIZE +  12*HOP_SIZE;
}
