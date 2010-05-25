/*
 * MATLAB Compiler: 4.10 (R2009a)
 * Date: Tue May 25 10:35:48 2010
 * Arguments: "-B" "macro_default" "-o" "PIVadvance4" "-W"
 * "WinMain:PIVadvance4" "-d"
 * "W:\matlab_projects\pivadvance\PIVadvance4\PIVadvance4\src" "-T" "link:exe"
 * "-v" "W:\matlab_projects\pivadvance\PIVadvance4\PIVadvance4.m" 
 */

#include "mclmcrrt.h"

#ifdef __cplusplus
extern "C" {
#endif
const unsigned char __MCC_PIVadvance4_session_key[] = {
    '5', 'B', '2', 'D', '1', '6', '2', '3', 'E', 'E', 'C', '9', 'E', '4', '7',
    '4', 'E', '0', '6', 'F', 'A', 'B', '2', '2', '2', 'A', '0', '3', '8', 'A',
    '0', 'C', '6', 'F', '2', '9', 'D', 'E', '8', '3', '7', 'D', 'F', '7', 'F',
    'E', '5', '6', '5', '0', 'D', '8', '0', 'E', 'D', 'D', 'A', '2', '0', 'A',
    '1', 'D', 'B', '0', '9', '4', '1', '2', '1', '9', '5', 'F', '5', 'E', '4',
    'A', '7', '8', 'F', 'B', '2', '3', 'F', 'B', '7', '9', 'B', '0', '6', 'B',
    '9', '8', 'D', '5', '3', '7', 'A', '1', 'D', 'B', '3', '6', '4', '4', '3',
    '0', 'E', '8', '3', '6', '8', '0', '2', '9', '3', '4', '5', 'C', '5', '5',
    'D', 'F', '8', '0', 'A', 'F', '6', '9', 'A', '5', 'D', 'F', '3', '3', '6',
    'C', '3', '4', '6', '8', '4', 'F', '8', 'F', '2', '6', '1', '3', '1', '2',
    'C', '3', 'E', '2', '3', '6', '5', '5', 'B', '3', '7', '8', '2', '3', 'B',
    'A', '5', 'F', '8', '7', 'B', '1', '8', '3', '7', '3', '0', '5', 'D', 'D',
    'E', '3', 'B', 'B', '1', '0', '3', '6', '1', 'B', 'F', '1', '0', '4', '7',
    'A', '0', '8', '9', '8', '8', '1', '6', '0', '0', '8', '0', '2', '2', '8',
    '3', '7', '1', '0', 'C', 'D', 'D', 'D', '5', '4', '7', 'C', '8', '2', '0',
    '5', '1', '2', '4', 'B', '9', 'E', 'F', '4', 'D', '3', '1', 'B', '9', '5',
    'E', 'C', 'A', '4', 'A', '7', 'B', '9', 'F', 'E', '1', 'C', '3', '9', '8',
    '4', '\0'};

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
    "toolbox/shared/dastudio/", "$TOOLBOXMATLABDIR/datamanager/",
    "toolbox/compiler/", "toolbox/images/colorspaces/",
    "toolbox/images/images/", "toolbox/images/iptformats/",
    "toolbox/images/iptutils/", "toolbox/shared/imageslib/",
    "toolbox/shared/spcuilib/", "toolbox/signal/signal/", "toolbox/stats/" };

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
  48,

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
  "PIVadvance4_342A98C48B0E2B308596140FF53D264B",

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


