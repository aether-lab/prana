/*
 * MATLAB Compiler: 4.10 (R2009a)
 * Date: Fri May 28 10:52:07 2010
 * Arguments: "-B" "macro_default" "-o" "PIVadvance4" "-W"
 * "WinMain:PIVadvance4" "-d"
 * "W:\matlab_projects\pivadvance\PIVadvance4\PIVadvance4\src" "-T" "link:exe"
 * "-v" "-N" "-p" "curvefit" "-p" "images" "-p" "signal" "-p" "distcomp"
 * "W:\matlab_projects\pivadvance\PIVadvance4\PIVadvance4.m" 
 */

#include "mclmcrrt.h"

#ifdef __cplusplus
extern "C" {
#endif
const unsigned char __MCC_PIVadvance4_session_key[] = {
    '0', '7', 'D', '5', '0', '0', '2', '7', '7', 'C', '0', '9', '1', '8', '3',
    '4', '7', '5', '6', '0', '8', '4', '0', '7', 'B', '7', '1', '0', 'B', '9',
    '2', '7', 'A', '0', 'B', '7', '7', '6', '8', '3', 'E', '6', '4', 'C', '7',
    '9', '9', 'C', '5', '4', 'D', 'C', '7', '8', '1', 'F', 'A', 'C', '5', '7',
    '2', '2', 'A', '9', '2', 'E', '5', 'D', 'E', 'E', 'F', 'C', '2', '0', 'F',
    '6', '0', '8', 'C', '4', 'E', 'F', '1', '2', '8', 'F', 'E', '8', 'E', '3',
    'C', '8', 'D', '2', '6', '1', 'A', '5', '4', 'D', '1', '6', 'C', '2', 'A',
    '7', '9', 'C', 'C', '6', '9', '4', 'A', '2', 'F', '7', 'E', '5', '0', '6',
    '3', '0', '9', '1', 'A', '1', 'E', '7', '0', 'D', '3', '9', '6', 'C', '7',
    '9', 'A', 'D', 'C', '7', 'B', 'A', '9', '0', 'D', '4', 'B', '8', 'A', '4',
    '3', 'E', 'D', 'B', '9', 'B', '7', '9', '1', '5', '5', '5', '5', 'F', '6',
    '8', 'D', '9', '0', '4', 'F', 'B', 'A', 'B', '1', 'C', '9', 'A', 'A', '9',
    'D', '9', '8', '1', '2', 'C', '3', '4', '5', '0', '7', '8', '6', '8', '3',
    'D', 'B', 'B', 'A', '3', '6', 'E', 'D', '4', '4', '2', '3', '3', '7', 'D',
    '0', '5', 'F', '6', '1', '7', '5', '7', 'D', '7', '1', '9', '0', '2', 'B',
    'A', '2', '1', '4', '3', '3', 'C', 'F', '0', '7', '7', '3', '3', '0', '0',
    'B', '0', '1', '4', 'E', '8', '9', 'A', '7', '0', '3', 'E', '1', 'B', '0',
    '7', '\0'};

const unsigned char __MCC_PIVadvance4_public_key[] = {
    '3', '0', '8', '1', '9', 'D', '3', '0', '0', 'D', '0', '6', '0', '9', '2',
    'A', '8', '6', '4', '8', '8', '6', 'F', '7', '0', 'D', '0', '1', '0', '1',
    '0', '1', '0', '5', '0', '0', '0', '3', '8', '1', '8', 'B', '0', '0', '3',
    '0', '8', '1', '8', '7', '0', '2', '8', '1', '8', '1', '0', '0', 'C', '4',
    '9', 'C', 'A', 'C', '3', '4', 'E', 'D', '1', '3', 'A', '5', '2', '0', '6',
    '5', '8', 'F', '6', 'F', '8', 'E', '0', '1', '3', '8', 'C', '4', '3', '1',
    '5', 'B', '4', '3', '1', '5', '2', '7', '7', 'E', 'D', '3', 'F', '7', 'D',
    'A', 'E', '5', '3', '0', '9', '9', 'D', 'B', '0', '8', 'E', 'E', '5', '8',
    '9', 'F', '8', '0', '4', 'D', '4', 'B', '9', '8', '1', '3', '2', '6', 'A',
    '5', '2', 'C', 'C', 'E', '4', '3', '8', '2', 'E', '9', 'F', '2', 'B', '4',
    'D', '0', '8', '5', 'E', 'B', '9', '5', '0', 'C', '7', 'A', 'B', '1', '2',
    'E', 'D', 'E', '2', 'D', '4', '1', '2', '9', '7', '8', '2', '0', 'E', '6',
    '3', '7', '7', 'A', '5', 'F', 'E', 'B', '5', '6', '8', '9', 'D', '4', 'E',
    '6', '0', '3', '2', 'F', '6', '0', 'C', '4', '3', '0', '7', '4', 'A', '0',
    '4', 'C', '2', '6', 'A', 'B', '7', '2', 'F', '5', '4', 'B', '5', '1', 'B',
    'B', '4', '6', '0', '5', '7', '8', '7', '8', '5', 'B', '1', '9', '9', '0',
    '1', '4', '3', '1', '4', 'A', '6', '5', 'F', '0', '9', '0', 'B', '6', '1',
    'F', 'C', '2', '0', '1', '6', '9', '4', '5', '3', 'B', '5', '8', 'F', 'C',
    '8', 'B', 'A', '4', '3', 'E', '6', '7', '7', '6', 'E', 'B', '7', 'E', 'C',
    'D', '3', '1', '7', '8', 'B', '5', '6', 'A', 'B', '0', 'F', 'A', '0', '6',
    'D', 'D', '6', '4', '9', '6', '7', 'C', 'B', '1', '4', '9', 'E', '5', '0',
    '2', '0', '1', '1', '1', '\0'};

static const char * MCC_PIVadvance4_matlabpath_data[] = 
  { "PIVadvance4/", "$TOOLBOXDEPLOYDIR/", "$TOOLBOXMATLABDIR/general/",
    "$TOOLBOXMATLABDIR/ops/", "$TOOLBOXMATLABDIR/lang/",
    "$TOOLBOXMATLABDIR/elmat/", "$TOOLBOXMATLABDIR/randfun/",
    "$TOOLBOXMATLABDIR/elfun/", "$TOOLBOXMATLABDIR/specfun/",
    "$TOOLBOXMATLABDIR/matfun/", "$TOOLBOXMATLABDIR/datafun/",
    "$TOOLBOXMATLABDIR/polyfun/", "$TOOLBOXMATLABDIR/funfun/",
    "$TOOLBOXMATLABDIR/sparfun/", "$TOOLBOXMATLABDIR/scribe/",
    "$TOOLBOXMATLABDIR/graph2d/", "$TOOLBOXMATLABDIR/graph3d/",
    "$TOOLBOXMATLABDIR/specgraph/", "$TOOLBOXMATLABDIR/graphics/",
    "$TOOLBOXMATLABDIR/uitools/", "$TOOLBOXMATLABDIR/strfun/",
    "$TOOLBOXMATLABDIR/imagesci/", "$TOOLBOXMATLABDIR/iofun/",
    "$TOOLBOXMATLABDIR/audiovideo/", "$TOOLBOXMATLABDIR/timefun/",
    "$TOOLBOXMATLABDIR/datatypes/", "$TOOLBOXMATLABDIR/verctrl/",
    "$TOOLBOXMATLABDIR/codetools/", "$TOOLBOXMATLABDIR/helptools/",
    "$TOOLBOXMATLABDIR/winfun/", "$TOOLBOXMATLABDIR/winfun/NET/",
    "$TOOLBOXMATLABDIR/demos/", "$TOOLBOXMATLABDIR/timeseries/",
    "$TOOLBOXMATLABDIR/hds/", "$TOOLBOXMATLABDIR/guide/",
    "$TOOLBOXMATLABDIR/plottools/", "toolbox/local/",
    "$TOOLBOXMATLABDIR/datamanager/", "toolbox/compiler/",
    "toolbox/images/colorspaces/", "toolbox/images/images/",
    "toolbox/images/iptformats/", "toolbox/images/iptutils/",
    "toolbox/shared/imageslib/", "toolbox/signal/signal/" };

static const char * MCC_PIVadvance4_classpath_data[] = 
  { "java/jar/toolbox/images.jar" };

static const char * MCC_PIVadvance4_libpath_data[] = 
  { "" };

static const char * MCC_PIVadvance4_app_opts_data[] = 
  { "" };

static const char * MCC_PIVadvance4_run_opts_data[] = 
  { "" };

static const char * MCC_PIVadvance4_warning_state_data[] = 
  { "off:MATLAB:dispatcher:nameConflict" };


mclComponentData __MCC_PIVadvance4_component_data = { 

  /* Public key data */
  __MCC_PIVadvance4_public_key,

  /* Component name */
  "PIVadvance4",

  /* Component Root */
  "",

  /* Application key data */
  __MCC_PIVadvance4_session_key,

  /* Component's MATLAB Path */
  MCC_PIVadvance4_matlabpath_data,

  /* Number of directories in the MATLAB Path */
  45,

  /* Component's Java class path */
  MCC_PIVadvance4_classpath_data,
  /* Number of directories in the Java class path */
  1,

  /* Component's load library path (for extra shared libraries) */
  MCC_PIVadvance4_libpath_data,
  /* Number of directories in the load library path */
  0,

  /* MCR instance-specific runtime options */
  MCC_PIVadvance4_app_opts_data,
  /* Number of MCR instance-specific runtime options */
  0,

  /* MCR global runtime options */
  MCC_PIVadvance4_run_opts_data,
  /* Number of MCR global runtime options */
  0,
  
  /* Component preferences directory */
  "PIVadvance4_A04E69EF73E225736BDB9E1DDBFF9CD5",

  /* MCR warning status data */
  MCC_PIVadvance4_warning_state_data,
  /* Number of MCR warning status modifiers */
  1,

  /* Path to component - evaluated at runtime */
  NULL

};

#ifdef __cplusplus
}
#endif


