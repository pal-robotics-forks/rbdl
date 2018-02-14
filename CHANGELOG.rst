^^^^^^^^^^^^^^^^^^^^^^^^^^
Changelog for package rbdl
^^^^^^^^^^^^^^^^^^^^^^^^^^

Forthcoming
-----------
* Merge branch 'as-urdfreader-linkskip' into 'erbium-devel'
  urdfreader: added functions to read URDF without given links
  See merge request control/rbdl!1
* Merge branch 'erbium-devel' into as-urdfreader-linkskip
* Dummy commit to trigger recursive testing
* Use quaternions instead of RPY angles to set joint transformations.
  Previous version seem to be introducing noise due to conversions back
  and forth (quaternion -> RPY -> rotation matrix).
* Model: added const getModelData().
* URDF reader: allow explicit specification of the kinamatic root link.
* URDF reader: minor changes in the logic
* Minor bugfix in URDF reader.
* Omit links in URDFModel: deleted corresponding functions from URDF reader
* urdfreader: code deduplication & cleanup
* Merge branch 'erbium-devel' into as-urdfreader-linkskip
  Conflicts:
  addons/urdfreader/urdfreader.cc
* urdfreader: drop unnecessary piece of code.
* Merge branch 'erbium-devel' into as-urdfreader-linkskip
  Conflicts:
  addons/urdfreader/urdfreader.cc
* Merge 'erbium-devel', cleanups, deduplications.
  Conflicts:
  addons/urdfreader/urdfreader.cc
  include/rbdl/addons/urdfreader/urdfreader.h
* Merge branch 'erbium-devel' into as-urdfreader-linkskip
  Conflicts:
  addons/urdfreader/urdfreader.cc
* urdfreader: added functions to read URDF without given links
  + some refactoring and partial formatting.
* Contributors: Hilario Tome, alexandersherikov

0.2.1 (2018-02-13)
------------------
* fixed compilation isnan
* Contributors: Hilario Tome

0.2.0 (2018-01-19)
------------------
* more templetization
* added rbdl parser function
* Merge branch 'erbium-devel' of gitlab:control/rbdl into erbium-devel
* fix template quaternion
* added extra parser
* changed rbdl root name for fixed floating base
* more bug fixes
* added proper root naming in fixed base rbdl
* fixed merge
* fixed critical bug in set body quaternion, the code was commented
* formating
* more templetization
* more templetization
* added specializations
* formating
* fixed utils
* more migration
* more migration
* fixed getter enum compile warking treated as error
* added better enum
* Merge branch 'dubnium-devel' into erbium-devel
* added get point angular acceleration and helper functions
* progres
* Merge branch 'dubnium-devel' into erbium-devel
* added helper util
* Added NO_TYPE floatingBaseType for grasping simulator
* unified utils
* More templetization
* Separated model data into a different header file
* Continue refactoring
* Broken commit, progress in having model as const
* Fixed bug
* Added model_data structure
* Templatized basic math operations
* Added coment
* Added utils
* Contributors: Adrià Roig, Hilario Tome, Hilario Tomé

0.1.1 (2016-10-14)
------------------
* Added conversion of mimic joints to fixed joints
* Contributors: Hilario Tome

0.1.0 (2016-10-05)
------------------
* Fixed cppecheck errors
* Update README.md
* Added gtests
* Updated to new RBDL version
* Update rbdl parser to parse urdf model
* Merge branch 'dubnium-devel' of gitlab:control/rbdl into dubnium-devel
* Added momentum computation
* Contributors: Hilario Tome

0.0.2 (2016-03-07)
------------------
* Added 2d floating base support
* Changed catkin package order in CMakeLists
* Aded various fixes and removed logging
* Removed unnecesary joint variable that created an allocation in update custom allocation
* Fixed allocation in joint
* Contributors: Hilario Tome

0.0.1 (2015-01-13)
------------------
* Release
