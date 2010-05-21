/*
 * MATLAB Compiler: 4.10 (R2009a)
 * Date: Fri May 21 11:08:54 2010
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
    '5', 'E', '6', '3', 'E', 'A', '6', 'C', 'E', '5', '2', 'A', 'A', '2', 'A',
    '0', 'A', 'A', 'A', '5', '7', '2', 'F', '5', '4', '8', 'B', '0', '8', 'E',
    'E', 'B', '5', 'E', '2', 'A', '9', 'C', 'F', '6', '2', '7', 'E', '1', '1',
    '7', '5', 'D', '1', 'D', 'B', '5', '3', '5', 'F', '8', 'D', 'B', 'F', '3',
    '8', '6', '1', 'E', '4', 'E', 'D', '0', 'F', '3', '4', '6', '3', 'B', 'C',
    '5', '8', '5', '8', '4', '4', '6', 'F', '4', 'D', 'C', '7', 'F', '3', '9',
    '7', '8', '7', '9', 'E', 'A', 'E', '9', '8', 'B', '1', '6', '7', '6', '9',
    '4', '9', '2', '1', 'F', 'B', '0', 'D', '9', '4', '1', '9', '8', '6', '7',
    '3', 'A', 'C', '1', '0', '0', '3', '1', '6', '3', '5', '6', 'B', '0', 'C',
    'D', 'E', 'D', '3', '9', '5', '2', '9', '6', '4', '7', '6', 'C', '7', '1',
    '7', '7', '5', 'D', 'C', 'A', 'B', 'C', '1', '0', '5', '6', '5', '3', '3',
    '5', 'F', 'C', '7', 'E', 'E', '0', 'B', '4', 'E', 'F', '7', 'C', 'B', '8',
    'C', '3', 'F', '7', 'B', '2', '2', 'E', '2', '1', '6', 'F', 'D', 'A', '9',
    '7', 'C', '8', '9', 'B', '1', 'A', 'D', 'E', '8', '0', 'F', '2', '4', '9',
    'F', '5', '7', '4', 'E', '7', 'C', 'F', '2', '5', 'D', 'D', 'E', '3', '8',
    '6', '4', '4', '5', '8', '6', 'A', 'D', '9', '9', '4', 'B', 'B', '6', '9',
    'A', '3', '2', '0', 'F', '0', 'B', '1', '2', '5', '8', 'B', 'F', 'F', '7',
    '1', '\0'};

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
    "$TOOLBOXMATLABDIR/winfun/", "$TOOLBOXMATLABDIR/winfun/net/",
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
  "PIVadvance4_C15D9D83987697A437A66793B0D7804B",

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


