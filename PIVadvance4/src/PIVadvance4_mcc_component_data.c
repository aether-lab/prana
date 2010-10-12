/*
 * MATLAB Compiler: 4.13 (R2010a)
 * Date: Tue Oct 12 14:30:25 2010
 * Arguments: "-B" "macro_default" "-o" "PIVadvance4" "-W"
 * "WinMain:PIVadvance4" "-T" "link:exe" "-d"
 * "W:\matlab_projects\pivadvance\PIVadvance4\PIVadvance4\src" "-w"
 * "enable:specified_file_mismatch" "-w" "enable:repeated_file" "-w"
 * "enable:switch_ignored" "-w" "enable:missing_lib_sentinel" "-w"
 * "enable:demo_license" "-v"
 * "W:\matlab_projects\pivadvance\PIVadvance4\PIVadvance4.m" "-a"
 * "W:\matlab_projects\pivadvance\PIVadvance4\AETHERlogo.mat" "-a"
 * "W:\matlab_projects\pivadvance\PIVadvance4\defaultsettings.mat" "-a"
 * "W:\matlab_projects\pivadvance\PIVadvance4\documentation" "-a"
 * "W:\matlab_projects\pivadvance\PIVadvance4\PIVadvance4.asv" "-a"
 * "W:\matlab_projects\pivadvance\PIVadvance4\PIVadvance4.fig" "-a"
 * "W:\matlab_projects\pivadvance\PIVadvance4\PIVadvance4.m~" "-a"
 * "W:\matlab_projects\pivadvance\PIVadvance4\PIVadvance4code.asv" "-a"
 * "W:\matlab_projects\pivadvance\PIVadvance4\PIVadvance4code.m" "-a"
 * "W:\matlab_projects\pivadvance\PIVadvance4\PIVadvance4code.m~" 
 */

#include "mclmcrrt.h"

#ifdef __cplusplus
extern "C" {
#endif
const unsigned char __MCC_PIVadvance4_session_key[] = {
    '6', 'C', '5', '6', '9', 'D', 'E', '9', 'F', 'D', '0', '7', '4', '9', 'D',
    'A', '1', 'C', '4', '6', 'D', '8', '4', 'E', '0', '3', '3', '7', 'A', 'A',
    '4', 'A', 'D', '2', '4', '9', '1', 'C', 'F', '9', 'D', 'F', '3', '1', '1',
    'A', '3', '9', '7', '1', '1', '0', '7', 'D', '2', '8', '2', 'C', '0', '5',
    '7', '6', 'E', 'F', '6', '4', '3', 'E', 'F', '5', 'C', '2', 'D', 'C', 'A',
    '9', '3', 'B', '8', '3', '1', 'F', '3', '4', 'E', 'A', '7', 'F', '0', '8',
    'B', 'E', 'D', 'B', '2', '2', '2', 'E', 'A', '9', '6', 'A', '8', '5', '2',
    'B', 'D', 'C', '3', '9', 'C', '7', '7', '1', 'C', '5', 'F', '4', 'D', 'A',
    '9', '3', '6', 'B', '6', 'C', 'B', '7', '9', 'C', '1', '6', '9', '9', '9',
    '8', 'E', 'F', '4', 'C', 'A', 'E', '9', '8', '4', 'A', '3', '5', 'E', 'F',
    'D', '4', 'B', '3', 'B', 'B', '1', '1', 'E', '9', '8', '5', '4', '5', 'F',
    '8', '0', '4', 'F', '5', '7', '2', 'E', '9', 'A', 'C', 'B', 'C', 'D', 'D',
    '9', '8', '3', 'B', 'B', 'B', 'A', 'A', '7', 'C', '7', '7', '1', 'B', '5',
    '1', 'A', '4', 'C', '3', '7', '5', '5', '2', 'E', '9', '6', '6', '5', '5',
    '1', '9', 'E', '2', '7', 'B', '1', '6', '9', 'F', 'D', '7', '9', '4', '6',
    '4', 'C', 'F', '6', 'D', '4', '3', '6', '4', '0', 'A', 'E', 'C', '1', '4',
    '8', '0', '0', 'B', '3', '6', '6', 'B', '4', 'A', 'D', 'A', '6', '6', '2',
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
  { "PIVadvance4/", "$TOOLBOXDEPLOYDIR/", "documentation/",
    "documentation/gettingstarted_files/", "$TOOLBOXMATLABDIR/general/",
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
    "toolbox/shared/optimlib/", "toolbox/distcomp/",
    "toolbox/distcomp/distcomp/", "toolbox/distcomp/mpi/",
    "toolbox/distcomp/parallel/", "toolbox/distcomp/parallel/lapack/",
    "toolbox/distcomp/parallel/util/", "toolbox/distcomp/lang/",
    "toolbox/shared/spcuilib/", "toolbox/images/colorspaces/",
    "toolbox/images/images/", "toolbox/images/imuitools/",
    "toolbox/images/iptformats/", "toolbox/images/iptutils/",
    "toolbox/shared/imageslib/", "toolbox/optim/optim/",
    "toolbox/signal/signal/", "toolbox/stats/" };

static const char * MCC_PIVadvance4_classpath_data[] = 
  { "java/jar/toolbox/images.jar" };

static const char * MCC_PIVadvance4_libpath_data[] = 
  { "bin/win32/" };

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
  59,

  /* Component's Java class path */
  MCC_PIVadvance4_classpath_data,
  /* Number of directories in the Java class path */
  1,

  /* Component's load library path (for extra shared libraries) */
  MCC_PIVadvance4_libpath_data,
  /* Number of directories in the load library path */
  1,

  /* MCR instance-specific runtime options */
  MCC_PIVadvance4_app_opts_data,
  /* Number of MCR instance-specific runtime options */
  0,

  /* MCR global runtime options */
  MCC_PIVadvance4_run_opts_data,
  /* Number of MCR global runtime options */
  0,
  
  /* Component preferences directory */
  "PIVadvance4_E703E3FF5275EB2A209C029388CD04AC",

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


