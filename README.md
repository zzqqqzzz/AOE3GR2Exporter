# AOE3GR2Exporter
Export AOE3 GR2 model and animation to FBX

## Requirements:
* Visual Studio 2019
* FBX SDK 2018.0
* Granny Common 2.11.8 SDK

## Usage:
* Extract gr2 files from AOE3 Art Assets
* `gr2fbxexporter.exe <modelfile> <outputfbxfile> [animfile1] [animfile2]`
* `gr2fbxexporter.exe <animfile1> <outputfbxfile>` works too.

### Example:
```
gr2fbxexporter.exe D:\age3w\Art\units\infantry\pikeman\pikeman_2age.gr2 pikeman_2age.fbx D:\age3w\Art\animation_library\infantry\Charge\pikemen_charge_attackB.gr2 D:\age3w\Art\animation_library\infantry\Volley\pikemen_volley_run.gr2
granny_file
-----------
File contains: 1 sections.
  Section 0: present (compressed)

Autodesk FBX SDK version 2018.0 Release (246225)
granny_file
-----------
File contains: 1 sections.
  Section 0: present (compressed)

granny_file
-----------
File contains: 1 sections.
  Section 0: present (compressed)

granny_file
-----------
File contains: 1 sections.
  Section 0: present (compressed)

FBX file format version 7.5.0

File pikeman_2age.fbx exported
```

Output sample in *Misc* folder

## Known issues:
Some gr2 files like rifleman.gr2 cannot be converted. Granny Viewer cannot open it either.

<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements
* [nwn2mdk by Arbos](https://github.com/Arbos/nwn2mdk)
* [FD Vertex Library by kangcliff](https://github.com/kangcliff/Age-of-Empires-III)
* [GrannyMeshDumper by Helia01](https://gyazo.com/1b325ce68d9293c94b2475ae62805304)
* [AoE3Ed by Ykkrosh](http://games.build-a.com/aoe3/files/)