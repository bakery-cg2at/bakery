Workflow
=========================================

The whole workflow of reverse mapping of molecules is divided into two parts as well as code. We will describe it on the example of dodecane.
The files are available in `example/dodecane` directory that is attached to the code source.

**Input files**

 - ``cg_conf.gro``: input CG coordinate file
 - ``dodecane_settings.xml``: input XML settings file
 - ``dodecane_single.gro``: input atomistic coordinate file

**Steps**

 1. ``cd example/dodecane``
 2. ``../../src/prepare_files.py --options dodecane_settings.xml`` - prepare hybrid configuration files
 3. ``../../src/start_backmapping.py --conf hyb_conf.gro --topology hyb_topology.top --alpha 0.00001`` - run reverse mapping procedure

In the end we will get two files:

  - ``sim0confout_final_aa_0.0001_12345.gro``
  - ``sim0topol_aa_0.0001_12345.top```.

