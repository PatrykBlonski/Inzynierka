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
 * @file panner_internal.h
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

#ifndef __PANNER_INTERNAL_H_INCLUDED__
#define __PANNER_INTERNAL_H_INCLUDED__

#include "panner.h"        /* Include header for this example */
#include "saf.h"           /* Main include header for SAF */
#include "saf_externals.h" /* To also include SAF dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#define FORCE_3D_LAYOUT /**< FLAG: Force 2D loudspeaker setups to also use 3D VBAP (i.e. with 2 virtual loudspeakers on the top/bottom) */
#if !defined(PANNER_FRAME_SIZE)
# if defined(FRAME_SIZE) /* Use the global framesize if it is specified: */
#  define PANNER_FRAME_SIZE ( FRAME_SIZE )          /**< Framesize, in time-domain samples */
# else /* Otherwise, the default framesize for this example is: */
#  define PANNER_FRAME_SIZE ( 480000 )                 /**< Framesize, in time-domain samples */
# endif
#endif
#define HOP_SIZE ( 128 )                            /**< STFT hop size */
#define HYBRID_BANDS ( HOP_SIZE + 5 )               /**< Number of frequency bands */
#define TIME_SLOTS ( PANNER_FRAME_SIZE / HOP_SIZE ) /**< Number of STFT timeslots */
#define MAX_COV_AVG_COEFF ( 0.45f )                   /**< Maximum supported covariance averaging coefficient  */
#define NUM_DISP_SLOTS ( 2 )                          /**< Number of display slots */


/* Checks: */
#if (PANNER_FRAME_SIZE % HOP_SIZE != 0)
# error "PANNER_FRAME_SIZE must be an integer multiple of HOP_SIZE"
#endif
    
/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/** Contains variables for scanning grids, and beamforming */
typedef struct _powermap_codecPars
{
    float* grid_dirs_deg;   /**< Spherical scanning grid directions, in degrees; FLAT: grid_nDirs x 2 */
    int grid_nDirs;         /**< Number of scanning directions */
    float* interp_dirs_deg; /**< 2D rectangular window interpolation directions, in degrees; FLAT: interp_nDirs x 2 */
    float* interp_table;    /**< Spherical->2D interpolation table; FLAT: interp_nDirs x grid_nDirs */
    int interp_nDirs;       /**< Number of interpolation directions */
    int interp_nTri;        /**< Number of triangles in the spherical triangulared grid */
    float* Y_grid[MAX_SH_ORDER];                 /**< real SH basis (real datatype); MAX_NUM_SH_SIGNALS x grid_nDirs */
    float_complex* Y_grid_cmplx[MAX_SH_ORDER];   /**< real SH basis (complex datatype); MAX_NUM_SH_SIGNALS x grid_nDirs */
    
}powermap_codecPars;

/**
 * Main structure for panner. Contains variables for audio buffers, afSTFT,
 * internal variables, flags, user parameters
 */
typedef struct _panner
{
    int new_nLoudpkrs;              /**< New number of loudspeakers in the array */
    int new_nSources;               /**< New number of inputs/sources */
    
    
    /* user parameters */
    int nSources;                   /**< Current number of inputs/sources */
    float src_dirs_deg[MAX_NUM_INPUTS][2]; /**< Current source directions */
    int nLoudpkrs;                  /**< Current number of loudspeakers in the array */
    float loudpkrs_dirs_deg[MAX_NUM_OUTPUTS][3]; /**< Current loudspeaker directions */




    /* FIFO buffers */
    int FIFO_idx;                   /**< FIFO buffer index */
    float inFIFO[MAX_NUM_SH_SIGNALS][PANNER_FRAME_SIZE]; /**< Input FIFO buffer */

    /* TFT */
    float** SHframeTD;              /**< time-domain SH input frame; #MAX_NUM_SH_SIGNALS x #POWERMAP_FRAME_SIZE */
    float_complex*** SHframeTF;     /**< time-frequency domain SH input frame; #HYBRID_BANDS x #MAX_NUM_SH_SIGNALS x #TIME_SLOTS */
    void* hSTFT;                    /**< afSTFT handle */
    void* spPWD;                    /**< sphPWD handle */
    float freqVector[HYBRID_BANDS]; /**< Frequency vector (filterbank centre frequencies) */
    float fs;                       /**< Host sample rate, in Hz*/

    /* internal */
    float_complex Cx[HYBRID_BANDS][MAX_NUM_SH_SIGNALS * MAX_NUM_SH_SIGNALS];     /**< covariance matrices per band */
    int new_masterOrder;            /**< New maximum/master SH analysis order (current value will be replaced by this after next re-init) */
    int dispWidth;                  /**< Number of pixels on the horizontal in the 2D interpolated powermap image */

    /* ana configuration */
    CODEC_STATUS codecStatus;       /**< see #CODEC_STATUS */
    PROC_STATUS procStatus;         /**< see #PROC_STATUS */
    float progressBar0_1;           /**< Current (re)initialisation progress, between [0..1] */
    char* progressBarText;          /**< Current (re)initialisation step, string */
    powermap_codecPars* pars;       /**< codec parameters */

    /* display */
    float* pmap;                    /**< grid_nDirs x 1 */
    float* prev_pmap;               /**< grid_nDirs x 1 */
    float* pmap_grid[NUM_DISP_SLOTS]; /**< powermap interpolated to grid; interp_nDirs x 1 */
    int dispSlotIdx;                /**< Current display slot */
    float pmap_grid_minVal;         /**< Current minimum value in pmap (used to normalise [0..1]) */
    float pmap_grid_maxVal;         /**< Current maximum value in pmap (used to normalise [0..1]) */
    int recalcPmap;                 /**< set this to 1 to generate a new powermap */
    int pmapReady;                  /**< 0: powermap not started yet, 1: powermap is ready for plotting*/

    /* User parameters */
    int masterOrder;                /**< Current maximum/master SH analysis order */
    int analysisOrderPerBand[HYBRID_BANDS]; /**< SH analysis order per frequency band */
    float pmapEQ[HYBRID_BANDS];     /**< Equalisation/weights per band */
    HFOV_OPTIONS HFOVoption;        /**< see #HFOV_OPTIONS */
    ASPECT_RATIO_OPTIONS aspectRatioOption; /**< see #ASPECT_RATIO_OPTIONS */
    float covAvgCoeff;              /**< Covariance matrix averaging coefficient, [0..1] */
    float pmapAvgCoeff;             /**< Powermap averaging coefficient, [0..1] */
    CH_ORDER chOrdering;            /**< Ambisonic channel order convention (see #CH_ORDER) */
    NORM_TYPES norm;                /**< Ambisonic normalisation convention (see #NORM_TYPES) */
    int maxInd;
    int loudNumber;
    float maxEnergy;
    
} panner_data;
     
typedef struct _beam {
    float theta;
    float phi;
    float sinTheta;
    float cosTheta;
    float sinPhi;
    float cosPhi;
    PROC_STATUS status;
} Direction;

#define PI 3.14159265358979323846f
#define THETA_STEPS 60
#define PHI_STEPS 45
/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/** Sets codec status (see #CODEC_STATUS enum) */
void panner_setCodecStatus(void* const hPan, CODEC_STATUS newStatus);
    
/**
 * Intialises the VBAP gain table used for panning.
 *
 * @note Call ambi_dec_initTFT() (if needed) before calling this function
 */
    
/**
 * Initialise the filterbank used by panner.
 *
 * @note Call this function before panner_initGainTables()
 */
void panner_initTFT(void* const hPan);
    
void panner_initAna(void* const hPan);


/**
 * Loads source directions from preset
 *
 * @param[in]  preset   See #SOURCE_CONFIG_PRESETS enum
 * @param[out] dirs_deg Source/loudspeaker directions
 * @param[out] newNCH   (&) new number of channels
 * @param[out] nDims    (&) estimate of the number of dimensions (2 or 3)
 */
void panner_loadSourcePreset(SOURCE_CONFIG_PRESETS preset,
                             float dirs_deg[MAX_NUM_INPUTS][2],
                             int* newNCH,
                             int* nDims);

/**
 * Loads source/loudspeaker directions from preset
 *
 * @param[in]  preset   See #LOUDSPEAKER_ARRAY_PRESETS enum
 * @param[out] dirs_deg Source/loudspeaker directions
 * @param[out] newNCH   (&) new number of channels
 * @param[out] nDims    (&) estimate of the number of dimensions (2 or 3)
 */



#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __PANNER_INTERNAL_H_INCLUDED__ */
