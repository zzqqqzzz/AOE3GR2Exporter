// gr2fbxexporter.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "granny.h"
#include "cgmath.h"
#include "gr.h"
#include <fbxsdk.h>
#ifdef FBXSDK_ENV_WIN
// On Windows platform need to include this to define  _msize()
#include <malloc.h>
#endif

using namespace std;

const double time_step = 1 / 30.0;

#ifdef IOS_REF
    #undef  IOS_REF
    #define IOS_REF (*(pManager->GetIOSettings()))
#endif

//fd code

typedef unsigned char fd_int8;
typedef unsigned short fd_int16;
typedef unsigned int fd_int32;
typedef float fd_float;

typedef struct fd_header fd_header;
typedef struct fd_packet_group fd_packet_group;
typedef struct fd_vertex_packet fd_vertex_packet;
typedef struct fd_vertex fd_vertex;
typedef struct fd_pwngt34332_vertex fd_pwngt34332_vertex;

#define FD_MAGIC 0xfdabcd01

struct fd_header {
    fd_int32 MagicValue;
    fd_int32 DataSize;
    fd_int32 CRC32;
    fd_int32 TotalVertexCount;
    fd_int32 GroupCount;
    fd_int32 GroupSize;
    fd_int32 PacketCount;
    fd_int32 PacketSize;
    fd_int32 VertexCount;
    fd_int32 VertexSize;
    fd_int32 PacketVertexCount; // sum of vertex count of all vertex packets, usually equals to VertexPacketCount * 4
    fd_int32 LargestBoneIndex; // largest bone index found in the vertices
    fd_int32 IncludeTangentVectors; // bool flag which specifies if the data originally contains tangent vectors
};

// Group of vertex packets weighed to the same bone
struct fd_packet_group {
    fd_int16 PacketCount;
    fd_int16 FirstVertexIndex;
    fd_int32 BoneIndex;
};

// Packet of 4 vertices
struct fd_vertex_packet {
    fd_float Position[3][4];
    fd_int32 Normal[4];
    fd_int32 Tangent[4];
    fd_float TextureCoordinates0[4][2];
};

// Vertex weighed to 2 bones
struct fd_vertex {
    fd_float Position[3];
    fd_int8 BoneWeights[2];
    fd_int8 BoneIndices[2];
    fd_int32 Normal;
    fd_int32 Tangent;
    fd_float TextureCoordinates0[2];
};

struct fd_pwngt34332_vertex {
    fd_float Position[3];
    fd_int8 BoneWeights[4];
    fd_int8 BoneIndices[4];
    fd_float Normal[3];
    fd_float Tangent[3];
    fd_float TextureCoordinates0[2];
};

#define FD_EXT_VEC(v, i) (fd_float)(*((const fd_int8*)&v + i)) / 63.5f - 1.0f

fd_int32 FDGetVertexCount(const void* Source) {
    return ((const fd_header*)Source)->TotalVertexCount;
}

fd_int32 FDGetVertices(const void* Source, granny_pwngt34332_vertex* DestVertices)
{
    const fd_int8* data = (const fd_int8*)Source;
    const fd_header* header = (const fd_header*)data;

    if (header->MagicValue != FD_MAGIC)
        return 0;
    //if (header->CRC32 != crc32(data + 12, header->DataSize))
    //	return 0;

    const fd_packet_group * groups = (const fd_packet_group*)(data + 64);
    const fd_vertex_packet * packets = (const fd_vertex_packet*)(data + 64 + header->GroupSize);
    const fd_vertex * vertices = (const fd_vertex*)(data + 64 + header->GroupSize + header->PacketSize);

    fd_int32 i = 0, p = 0;
    for (fd_int32 g = 0; g < header->GroupCount; ++g) {
        for (fd_int32 h = 0; h < groups[g].PacketCount; ++h) {
            for (fd_int32 v = 0; v < 4; ++v) {
                DestVertices[i].Position[0] = packets[p].Position[0][v];
                DestVertices[i].Position[1] = packets[p].Position[1][v];
                DestVertices[i].Position[2] = packets[p].Position[2][v];
                DestVertices[i].BoneWeights[0] = 255;
                DestVertices[i].BoneIndices[0] = groups[g].BoneIndex;
                DestVertices[i].BoneWeights[1] = 0;
                DestVertices[i].BoneIndices[1] = 0;
                DestVertices[i].BoneWeights[2] = 0;
                DestVertices[i].BoneIndices[2] = 0;
                DestVertices[i].BoneWeights[3] = 0;
                DestVertices[i].BoneIndices[3] = 0;
                DestVertices[i].Normal[0] = FD_EXT_VEC(packets[p].Normal[v], 0);
                DestVertices[i].Normal[1] = FD_EXT_VEC(packets[p].Normal[v], 1);
                DestVertices[i].Normal[2] = FD_EXT_VEC(packets[p].Normal[v], 2);
                DestVertices[i].Tangent[0] = FD_EXT_VEC(packets[p].Tangent[v], 0);
                DestVertices[i].Tangent[1] = FD_EXT_VEC(packets[p].Tangent[v], 1);
                DestVertices[i].Tangent[2] = FD_EXT_VEC(packets[p].Tangent[v], 2);
                DestVertices[i].UV[0] = packets[p].TextureCoordinates0[v][0];
                DestVertices[i].UV[1] = packets[p].TextureCoordinates0[v][1];

                ++i;
            }
            ++p;
        }

    }

    for (fd_int32 v = 0; v < header->VertexCount; ++v) {
        DestVertices[i].Position[0] = vertices[v].Position[0];
        DestVertices[i].Position[1] = vertices[v].Position[1];
        DestVertices[i].Position[2] = vertices[v].Position[2];
        DestVertices[i].BoneWeights[0] = vertices[v].BoneWeights[0];
        DestVertices[i].BoneIndices[0] = vertices[v].BoneIndices[0];
        DestVertices[i].BoneWeights[1] = vertices[v].BoneWeights[1];
        DestVertices[i].BoneIndices[1] = vertices[v].BoneIndices[1];
        DestVertices[i].BoneWeights[2] = 0;
        DestVertices[i].BoneIndices[2] = 0;
        DestVertices[i].BoneWeights[3] = 0;
        DestVertices[i].BoneIndices[3] = 0;
        DestVertices[i].Normal[0] = FD_EXT_VEC(vertices[v].Normal, 0);
        DestVertices[i].Normal[1] = FD_EXT_VEC(vertices[v].Normal, 1);
        DestVertices[i].Normal[2] = FD_EXT_VEC(vertices[v].Normal, 2);
        DestVertices[i].Tangent[0] = FD_EXT_VEC(vertices[v].Tangent, 0);
        DestVertices[i].Tangent[1] = FD_EXT_VEC(vertices[v].Tangent, 1);
        DestVertices[i].Tangent[2] = FD_EXT_VEC(vertices[v].Tangent, 2);
        DestVertices[i].UV[0] = vertices[v].TextureCoordinates0[0];
        DestVertices[i].UV[1] = vertices[v].TextureCoordinates0[1];

        ++i;
    }

    return 1;
}

// end fd code

void PrintFileStatistics(granny_file*);

void InitializeSdkObjects(FbxManager*& pManager, FbxScene*& pScene)
{
    //The first thing to do is to create the FBX Manager which is the object allocator for almost all the classes in the SDK
    pManager = FbxManager::Create();
    if (!pManager)
    {
        FBXSDK_printf("Error: Unable to create FBX Manager!\n");
        exit(1);
    }
    else FBXSDK_printf("Autodesk FBX SDK version %s\n", pManager->GetVersion());

    //Create an IOSettings object. This object holds all import/export settings.
    FbxIOSettings* ios = FbxIOSettings::Create(pManager, IOSROOT);
    pManager->SetIOSettings(ios);

    //Load plugins from the executable directory (optional)
    FbxString lPath = FbxGetApplicationDirectory();
    pManager->LoadPluginsDirectory(lPath.Buffer());

    //Create an FBX scene. This object holds most objects imported/exported from/to files.
    pScene = FbxScene::Create(pManager, "My Scene");
    if (!pScene)
    {
        FBXSDK_printf("Error: Unable to create FBX scene!\n");
        exit(1);
    }
}

static FbxVector4 quat_to_euler(FbxQuaternion& q)
{
    FbxAMatrix m;
    m.SetQ(q);
    return m.GetR();
}

static FbxNode* create_node(FbxScene* scene, FbxMesh* mesh, const char* name)
{
    auto node = FbxNode::Create(scene, name);
    node->SetNodeAttribute(mesh);
    // Blender default import settings rotate 90 degrees around x-axis and
    // scale by 0.01. We appply the inverse transform to the FBX node.
    node->LclRotation.Set(FbxDouble3(0, 0, 0));
    node->LclScaling.Set(FbxDouble3(1, 1, 1));

    return node;
}

std::vector<FbxCluster*> create_clusters(granny_model* Model, FbxNode** skeletons, FbxScene* pScene, vector<granny_pwn343_vertex> Vertices, int VertexCount, granny_mesh* Mesh)
{
    granny_mesh_binding* meshbinding = GrannyNewMeshBinding(Mesh, Model->Skeleton, Model->Skeleton);
    const int* FromBoneIndices = GrannyGetMeshBindingFromBoneIndices(meshbinding);
    int const NumMeshBones = GrannyGetMeshBindingBoneCount(meshbinding);

    std::vector<FbxCluster*> clusters(Model->Skeleton->BoneCount);
    for (unsigned i = 0; i < clusters.size(); ++i) {
        auto fbx_bone = skeletons[i];
        clusters[i] = FbxCluster::Create(pScene, "");
        clusters[i]->SetLink(fbx_bone);
        clusters[i]->SetLinkMode(FbxCluster::eNormalize);
        FbxAMatrix m;
        m.SetT(FbxVector4(0, 0, 0));
        m.SetR(FbxVector4(0, 0, 0));
        m.SetS(FbxVector4(1, 1, 1));
        clusters[i]->SetTransformMatrix(m);
        clusters[i]->SetTransformLinkMatrix(fbx_bone->EvaluateGlobalTransform());
    }

    for (unsigned i = 0; i < VertexCount; ++i) {
        granny_pwn343_vertex& Vert = Vertices[i];
        
        for (unsigned j = 0; j < 4; ++j) {
            if (Vert.BoneWeights[j] > 0) {
                auto bone_index = FromBoneIndices[Vert.BoneIndices[j]];
                if (bone_index < clusters.size())
                    clusters[bone_index]->AddControlPointIndex(i, Vert.BoneWeights[j]);
                else
                    cout << "  WARNING: bone index out of bounds (" << unsigned(bone_index) << " >= " << clusters.size() << ")\n";
            }
        }
    }

    GrannyFreeMeshBinding(meshbinding);
    return clusters;
}

bool CreateScene(FbxScene* pScene, granny_file_info* FileInfo) {
    if (FileInfo->ModelCount == 0) {
        return true;
    }
    granny_model* Model = FileInfo->Models[0];
    // set bone
    int bonecount = Model->Skeleton->BoneCount;
    FbxNode** skeletons = (FbxNode * *)malloc(sizeof(FbxNode *) * bonecount);;
    for (int i = 0; i < Model->Skeleton->BoneCount; i++) {
        granny_bone bone = Model->Skeleton->Bones[i];
        FbxNode* pNode = FbxNode::Create(pScene, bone.Name);
        FbxSkeleton* lSkeleton = FbxSkeleton::Create(pScene, "");
        if (bone.ParentIndex == -1) {
            lSkeleton->SetSkeletonType(FbxSkeleton::eRoot);
        }
        else {
            lSkeleton->SetSkeletonType(FbxSkeleton::eLimbNode);
        }

        if (bone.LocalTransform.Flags & 0x1) {
            pNode->LclTranslation.Set(FbxDouble3(bone.LocalTransform.Position[0], bone.LocalTransform.Position[1], bone.LocalTransform.Position[2]));
        }
        else {
            pNode->LclTranslation.Set(FbxDouble3(0, 0, 0));
        }
        if (bone.LocalTransform.Flags & 0x2) {
            FbxQuaternion rotation(bone.LocalTransform.Orientation[0],
                bone.LocalTransform.Orientation[1],
                bone.LocalTransform.Orientation[2],
                bone.LocalTransform.Orientation[3]);
            pNode->LclRotation.Set(quat_to_euler(rotation));
        }
        else {
            pNode->LclRotation.Set(FbxDouble3(0, 0, 0));
        }
        pNode->LclScaling.Set(FbxDouble3(1,1,1));
        pNode->SetNodeAttribute(lSkeleton);
        pNode->TranslationActive = true;
        pNode->TranslationMinX = true;
        pNode->TranslationMinY = true;
        pNode->TranslationMinZ = true;

        skeletons[i] = pNode;
        if (bone.ParentIndex == -1) {
            pScene->GetRootNode()->AddChild(pNode);
        }
        else {
            skeletons[bone.ParentIndex]->AddChild(pNode);
        }
    }

    // set mesh
    for (granny_int32x i = 0; i < FileInfo->MeshCount; i++)
    {
        granny_mesh* Mesh = FileInfo->Meshes[i];
        if (!Mesh)
            continue;

        granny_variant Ignore;

        // convert FD to granny
        if (GrannyFindMatchingMember(GrannyGetMeshVertexType(Mesh),
            GrannyGetMeshVertices(Mesh),
            "FD",
            &Ignore))
        {
            granny_pwngt34332_vertex* vertices;
            int vertices_count = FDGetVertexCount(Mesh->PrimaryVertexData->Vertices);
            vertices = new granny_pwngt34332_vertex[vertices_count];

            if (!FDGetVertices(Mesh->PrimaryVertexData->Vertices, vertices))
                return false;

            Mesh->PrimaryVertexData->Vertices = (granny_uint8*)vertices;
            Mesh->PrimaryVertexData->VertexType = GrannyPWNGT34332VertexType;
            Mesh->PrimaryVertexData->VertexCount = vertices_count;
        }

        std::vector<int> Positions;
        
        if (GrannyFindMatchingMember(GrannyGetMeshVertexType(Mesh),
            GrannyGetMeshVertices(Mesh),
            GrannyVertexPositionName,
            &Ignore))
        {
            granny_data_type_definition PosType[] = {
                { GrannyReal32Member, GrannyVertexPositionName, 0, 3 },
                { GrannyEndMember }
            };
            Positions.resize(GrannyGetMeshVertexCount(Mesh) * 3);
            GrannyCopyMeshVertices(Mesh, PosType, &Positions[0]);
        }


        std::vector<float> Normals;

        if (GrannyFindMatchingMember(GrannyGetMeshVertexType(Mesh),
            GrannyGetMeshVertices(Mesh),
            GrannyVertexNormalName,
            &Ignore))
        {
            granny_data_type_definition NormalType[] = {
                { GrannyReal32Member, GrannyVertexNormalName, 0, 3 },
                { GrannyEndMember }
            };
            Normals.resize(GrannyGetMeshVertexCount(Mesh) * 3);
            GrannyCopyMeshVertices(Mesh, NormalType, &Normals[0]);
        }

        granny_int32 VertexBufferSize = GrannyGetMeshVertexCount(Mesh) * sizeof(granny_pnt332_vertex);
        void* VertexBufferPointer = malloc(VertexBufferSize);
        GrannyCopyMeshVertices(Mesh, GrannyPNT332VertexType, VertexBufferPointer);
        
        FbxMesh* mesh = FbxMesh::Create(pScene, Mesh->Name);
        int positionCount = GrannyGetMeshVertexCount(Mesh);
        mesh->InitControlPoints(positionCount);
        for (size_t p = 0; p < positionCount; p++) {
            granny_real32 p1 = ((granny_pnt332_vertex*)VertexBufferPointer)[p].Position[0];
            float* pp1 = &p1;
            granny_real32 p2 = ((granny_pnt332_vertex*)VertexBufferPointer)[p].Position[1];
            float* pp2 = &p2;
            FbxVector4 v(((granny_pnt332_vertex*)VertexBufferPointer)[p].Position[0], ((granny_pnt332_vertex*)VertexBufferPointer)[p].Position[1], ((granny_pnt332_vertex*)VertexBufferPointer)[p].Position[2]);
            
            mesh->GetControlPoints()[p] = v;
        }

        FbxGeometryElementNormal* lElementNormal = mesh->CreateElementNormal();
        lElementNormal->SetMappingMode(FbxGeometryElement::eByControlPoint);
        // Set the normal values for every control point.

        lElementNormal->SetReferenceMode(FbxGeometryElement::eDirect);

        for (size_t n = 0; n < positionCount; n++) {
            FbxVector4 lNormalZPos(((granny_pnt332_vertex*)VertexBufferPointer)[n].Normal[0], ((granny_pnt332_vertex*)VertexBufferPointer)[n].Normal[1], ((granny_pnt332_vertex*)VertexBufferPointer)[n].Normal[2]);
            lElementNormal->GetDirectArray().Add(lNormalZPos);
        }

        FbxGeometryElementUV* element_uv = mesh->CreateElementUV("UVMap");
        // One UV for each vertex.
        element_uv->SetMappingMode(FbxGeometryElement::eByControlPoint);
        // Map information in DirectArray.
        element_uv->SetReferenceMode(FbxGeometryElement::eDirect);

        for (uint32_t n = 0; n < positionCount; n++) {
            FbxVector2 uv((((granny_pnt332_vertex*)VertexBufferPointer)[n].UV[0]), (1 - ((granny_pnt332_vertex*)VertexBufferPointer)[n].UV[1]));
            element_uv->GetDirectArray().Add(uv);
        }

        delete(VertexBufferPointer);

        granny_int32 TriangleCount = GrannyGetMeshTriangleCount(Mesh);
        granny_int32 BytesPerIndex = GrannyGetMeshBytesPerIndex(Mesh);

        int biggestIndex = 0;
        void* Indices = GrannyGetMeshIndices(Mesh);
        for (int k = 0; k < TriangleCount; k++) {
            mesh->BeginPolygon(-1, -1, -1, false);

            if (BytesPerIndex == 4){
                mesh->AddPolygon(((int*)Indices)[k * 3]);
                mesh->AddPolygon(((int*)Indices)[k * 3 + 1]);
                mesh->AddPolygon(((int*)Indices)[k * 3 + 2]);
                mesh->EndPolygon();
            }
            else if (BytesPerIndex == 2) {
                if (biggestIndex < ((short*)Indices)[k * 3]) {
                    biggestIndex = ((short*)Indices)[k * 3];
                }
                if (biggestIndex < ((short*)Indices)[k * 3+1]) {
                    biggestIndex = ((short*)Indices)[k * 3+1];
                }
                if (biggestIndex < ((short*)Indices)[k * 3+2]) {
                    biggestIndex = ((short*)Indices)[k * 3+2];
                }
                mesh->AddPolygon(((short*)Indices)[k * 3]);
                mesh->AddPolygon(((short*)Indices)[k * 3 + 1]);
                mesh->AddPolygon(((short*)Indices)[k * 3 + 2]);
                mesh->EndPolygon();

            }
        }

        FbxNode* node = create_node(pScene, mesh, Mesh->Name);
        pScene->GetRootNode()->AddChild(node);


        int VertexCount = GrannyGetMeshVertexCount(Mesh);

        vector<granny_pwn343_vertex> Vertices(VertexCount);
        GrannyCopyMeshVertices(Mesh, GrannyPWN343VertexType, &Vertices[0]);
        GrannyOneNormalizeWeights(VertexCount, GrannyPWN343VertexType, &Vertices[0]);

        // create skin
        FbxSkin* fbx_skin = FbxSkin::Create(pScene, "");
        fbx_skin->SetSkinningType(FbxSkin::eRigid);
        auto clusters = create_clusters(Model, skeletons, pScene, Vertices, VertexCount, Mesh);
        for (unsigned i = 0; i < clusters.size(); ++i)
            fbx_skin->AddCluster(clusters[i]);

        mesh->AddDeformer(fbx_skin);

    }
    return true;
}

bool node_is_animated(FbxNode* node, granny_transform_track& transform_track)
{
    return strcmp(node->GetName(), transform_track.Name) == 0 &&
        (node->GetChildCount() == 0 || strcmp(node->GetName(), node->GetChild(0)->GetName()) != 0);
}

std::vector<float> padded_knots(const std::vector<float>& knots, unsigned degree)
{
    std::vector<float> v;

    v.push_back(0);

    for (auto t : knots)
        v.push_back(t);

    for (unsigned i = 1; i <= degree; ++i) {
        v[i] = 0;
        v.push_back(knots.back());
    }

    return v;
}

// De Boor's algorithm to evaluate a B-spline
Vector3<float> de_boor_position(unsigned k, const std::vector<float>& knots, const std::vector<Vector4<float>>& controls, float t)
{
    unsigned i = k;

    while (i < knots.size() - k - 1 && knots[i] <= t)
        ++i;

    i = i - 1;

    std::vector<Vector3<float>> d;

    for (unsigned j = 0; j <= k; ++j) {
        auto r = controls[j + i - k];
        d.emplace_back(r.x, r.y, r.z);
    }

    for (unsigned r = 1; r <= k; ++r) {
        for (unsigned j = k; j >= r; --j) {
            float alpha = 0;
            if (knots[j + 1 + i - r] != knots[j + i - k])
                alpha = (t - knots[j + i - k]) / (knots[j + 1 + i - r] - knots[j + i - k]);
            d[j].x = d[j - 1].x * (1 - alpha) + d[j].x * alpha;
            d[j].y = d[j - 1].y * (1 - alpha) + d[j].y * alpha;
            d[j].z = d[j - 1].z * (1 - alpha) + d[j].z * alpha;
        }
    }

    return d[k];
}

Vector3<float> de_boor_position(unsigned degree, float t, GR2_curve_view & v)
{
    if (v.knots().size() == 1 || v.knots().size() < degree + 1)
        return Vector3<float>(v.controls()[0].x, v.controls()[0].y, v.controls()[0].z);

    return de_boor_position(degree, padded_knots(v.knots(), degree), v.controls(), t);
}

// De Boor's algorithm to evaluate a B-spline
FbxQuaternion de_boor_rotation(unsigned k, const std::vector<float> & knots, const std::vector<Vector4<float>> & controls, float t)
{
    unsigned i = k;

    while (i < knots.size() - k - 1 && knots[i] <= t)
        ++i;

    i = i - 1;

    std::vector<FbxQuaternion> d;

    for (unsigned j = 0; j <= k; ++j) {
        auto r = controls[j + i - k];
        d.emplace_back(r.x, r.y, r.z, r.w);
    }

    for (unsigned r = 1; r <= k; ++r) {
        for (unsigned j = k; j >= r; --j) {
            float alpha = 1;
            if (knots[j + 1 + i - r] != knots[j + i - k])
                alpha = (t - knots[j + i - k]) / (knots[j + 1 + i - r] - knots[j + i - k]);
            d[j] = d[j - 1].Slerp(d[j], alpha);
        }
    }

    return d[k];
}

FbxQuaternion de_boor_rotation(unsigned degree, float t, GR2_curve_view & v)
{
    if (v.knots().size() == 1 || v.knots().size() < degree + 1)
        return FbxQuaternion(v.controls()[0].x, v.controls()[0].y, v.controls()[0].z, v.controls()[0].w);

    return de_boor_rotation(degree, padded_knots(v.knots(), degree), v.controls(), t);
}

void add_position_keyframe(FbxAnimCurve * curvex, FbxAnimCurve * curvey,
    FbxAnimCurve * curvez, GR2_curve_view & view, float t)
{
    auto p = de_boor_position(view.degree(), t, view);

    FbxTime time;
    time.SetSecondDouble(t);

    auto k = curvex->KeyAdd(time);
    curvex->KeySetInterpolation(k, FbxAnimCurveDef::eInterpolationLinear);
    curvex->KeySetValue(k, p.x);

    k = curvey->KeyAdd(time);
    curvey->KeySetInterpolation(k, FbxAnimCurveDef::eInterpolationLinear);
    // removed / 100
    curvey->KeySetValue(k, p.y);

    k = curvez->KeyAdd(time);
    curvez->KeySetInterpolation(k, FbxAnimCurveDef::eInterpolationLinear);
    curvez->KeySetValue(k, p.z);
    curvez->KeySetValue(k, p.z);
}

void create_anim_position(FbxNode * node, FbxAnimLayer * anim_layer, granny_animation * anim,
    granny_transform_track& transform_track)
{
    GR2_curve_view view(transform_track.PositionCurve);

    if (view.knots().empty())
        return;

    auto curvex = node->LclTranslation.GetCurve(anim_layer, FBXSDK_CURVENODE_COMPONENT_X, true);
    auto curvey = node->LclTranslation.GetCurve(anim_layer, FBXSDK_CURVENODE_COMPONENT_Y, true);
    auto curvez = node->LclTranslation.GetCurve(anim_layer, FBXSDK_CURVENODE_COMPONENT_Z, true);

    curvex->KeyModifyBegin();
    curvey->KeyModifyBegin();
    curvez->KeyModifyBegin();

    for (double i = 0, t = 0; t < anim->Duration + time_step / 2; ++i, t = i * time_step)
        add_position_keyframe(curvex, curvey, curvez, view, float(t));

    curvex->KeyModifyEnd();
    curvey->KeyModifyEnd();
    curvez->KeyModifyEnd();
}

void add_rotation_keyframe(FbxAnimCurve * curvex, FbxAnimCurve * curvey,
    FbxAnimCurve * curvez, GR2_curve_view & view, float t)
{
    auto r = de_boor_rotation(view.degree(), t, view);
    auto e = quat_to_euler(r);

    FbxTime time;
    time.SetSecondDouble(t);

    auto k = curvex->KeyAdd(time);
    curvex->KeySetInterpolation(k, FbxAnimCurveDef::eInterpolationLinear);
    curvex->KeySetValue(k, float(e[0]));

    k = curvey->KeyAdd(time);
    curvey->KeySetInterpolation(k, FbxAnimCurveDef::eInterpolationLinear);
    curvey->KeySetValue(k, float(e[1]));

    k = curvez->KeyAdd(time);
    curvez->KeySetInterpolation(k, FbxAnimCurveDef::eInterpolationLinear);
    curvez->KeySetValue(k, float(e[2]));
}

void create_anim_rotation(FbxNode * node, FbxAnimLayer * anim_layer, granny_animation * anim, granny_transform_track& transform_track)
{
    GR2_curve_view view(transform_track.OrientationCurve);

    if (view.knots().empty())
        return;

    auto curvex = node->LclRotation.GetCurve(anim_layer, FBXSDK_CURVENODE_COMPONENT_X, true);
    auto curvey = node->LclRotation.GetCurve(anim_layer, FBXSDK_CURVENODE_COMPONENT_Y, true);
    auto curvez = node->LclRotation.GetCurve(anim_layer, FBXSDK_CURVENODE_COMPONENT_Z, true);

    curvex->KeyModifyBegin();
    curvey->KeyModifyBegin();
    curvez->KeyModifyBegin();

    for (double i = 0, t = 0; t < anim->Duration + time_step / 2; ++i, t = i * time_step)
        add_rotation_keyframe(curvex, curvey, curvez, view, float(t));

    curvex->KeyModifyEnd();
    curvey->KeyModifyEnd();
    curvez->KeyModifyEnd();
}

std::pair<std::vector<float>, std::vector<float>> scaleshear_curve_view(granny_transform_track& transform_track)
{
    std::vector<float> knots;
    std::vector<float> controls;
    granny_curve_data_header* header = (granny_curve_data_header*)transform_track.ScaleShearCurve.CurveData.Object;
    if (header->Format == DaConstant32f) {
        knots.push_back(0);
        auto data = (granny_curve_data_da_constant32f*)transform_track.ScaleShearCurve.CurveData.Object;
        for (int i = 0; i < 9; ++i)
            controls.push_back(data->Controls[i]);
    }
    else if (header->Format == DaK16uC16u) {
        auto data = (granny_curve_data_da_k16u_c16u*)transform_track.ScaleShearCurve.CurveData.Object;
        GR2_DaK16uC16u_view view(*data);
        return { view.knots(), view.controls() };
    }
    else if (header->Format == DaK32fC32f) {
        auto data = (granny_curve_data_da_k32f_c32f*)transform_track.ScaleShearCurve.CurveData.Object;
        for (int i = 0; i < data->KnotCount; ++i)
            knots.push_back(data->Knots[i]);
        for (int i = 0; i < data->KnotCount; ++i)
            controls.push_back(data->Controls[i]);
    }

    return { knots, controls };
}

void add_scaling_keyframe(FbxAnimCurve * curvex, FbxAnimCurve * curvey,
    FbxAnimCurve * curvez, std::vector<float> & knots, std::vector<float> & controls,
    unsigned i)
{
    FbxTime time;
    time.SetSecondDouble(knots[i]);

    auto k = curvex->KeyAdd(time);
    curvex->KeySetInterpolation(k, FbxAnimCurveDef::eInterpolationLinear);
    curvex->KeySetValue(k, controls[i * 9 + 0]);

    k = curvey->KeyAdd(time);
    curvey->KeySetInterpolation(k, FbxAnimCurveDef::eInterpolationLinear);
    curvey->KeySetValue(k, controls[i * 9 + 4]);

    k = curvez->KeyAdd(time);
    curvez->KeySetInterpolation(k, FbxAnimCurveDef::eInterpolationLinear);
    curvez->KeySetValue(k, controls[i * 9 + 8]);
}

void create_anim_scaling(FbxNode * node, FbxAnimLayer * anim_layer, granny_animation* anim, granny_transform_track& transform_track)
{
    std::pair<std::vector<float>, std::vector<float>> pair = scaleshear_curve_view(transform_track);
    std::vector<float> knots = pair.first;
    std::vector<float> controls = pair.second;

    if (knots.empty())
        return;

    auto curvex = node->LclScaling.GetCurve(anim_layer, FBXSDK_CURVENODE_COMPONENT_X, true);
    auto curvey = node->LclScaling.GetCurve(anim_layer, FBXSDK_CURVENODE_COMPONENT_Y, true);
    auto curvez = node->LclScaling.GetCurve(anim_layer, FBXSDK_CURVENODE_COMPONENT_Z, true);

    curvex->KeyModifyBegin();
    curvey->KeyModifyBegin();
    curvez->KeyModifyBegin();

    for (unsigned i = 0; i < knots.size(); ++i)
        add_scaling_keyframe(curvex, curvey, curvez, knots, controls, i);

    curvex->KeyModifyEnd();
    curvey->KeyModifyEnd();
    curvez->KeyModifyEnd();
}

void create_animation(FbxNode* node, FbxAnimLayer* anim_layer,
    granny_animation* anim, granny_transform_track& transform_track)
{
    create_anim_position(node, anim_layer, anim, transform_track);
    create_anim_rotation(node, anim_layer, anim, transform_track);
    create_anim_scaling(node, anim_layer, anim, transform_track);
}

FbxNode* create_animation_pivot(FbxNode* node, granny_track_group* track_group)
{
    string pivot_name = track_group->Name;
    if (strcmp(node->GetName(), track_group->Name) == 0)
        pivot_name += ".PIVOT";

    auto pivot = node->GetScene()->FindNodeByName(pivot_name.c_str());

    if (!pivot) {
        pivot = FbxNode::Create(node->GetScene(), pivot_name.c_str());
        pivot->LclRotation.Set(FbxDouble3(0, 0, 0));
        pivot->LclScaling.Set(FbxDouble3(1, 1, 1));
        node->GetScene()->GetRootNode()->AddChild(pivot);
    }

    return pivot;
}

void parent_to_animation_pivot(FbxNode* node, granny_track_group* track_group)
{
    if (!node->GetMesh())
        return; // Animation pivot is only for meshes

    auto pivot = create_animation_pivot(node, track_group);

    // Make node child of pivot
    node->GetParent()->RemoveChild(node);
    pivot->AddChild(node);

    // Reset rotation and scaling since the node is not root anymore
    node->LclRotation.Set(FbxDouble3(0, 0, 0));
    node->LclScaling.Set(FbxDouble3(1, 1, 1));
}

bool export_animation(FbxNode* node, FbxAnimLayer* anim_layer,
    granny_animation* anim, granny_track_group* track_group, granny_transform_track& transform_track)
{
    if (node_is_animated(node, transform_track)) {
        create_animation(node, anim_layer, anim, transform_track);
        parent_to_animation_pivot(node, track_group);
        return true;
    }

    for (int i = 0; i < node->GetChildCount(); ++i)
        if (export_animation(node->GetChild(i), anim_layer, anim, track_group, transform_track))
            return true;

    return false;
}

void export_animation(FbxScene* scene, granny_animation* anim, granny_track_group* track_group,
    FbxAnimLayer* anim_layer)
{
    for (int i = 0; i < track_group->TransformTrackCount; ++i) {
        export_animation(scene->GetRootNode(), anim_layer, anim, track_group, track_group->TransformTracks[i]);
    }
}

bool AddAnim(FbxScene* pScene, const char* filename) {
    granny_file* GrannyAnimFile = GrannyReadEntireFile(filename);
    if (GrannyAnimFile == 0)
    {
        printf("Error: unable to load %s as a Granny file.\n", filename);
        return EXIT_FAILURE;
    }
    else
    {
        PrintFileStatistics(GrannyAnimFile);
    }

    granny_file_info* AnimFileInfo = GrannyGetFileInfo(GrannyAnimFile);
    if (AnimFileInfo == NULL)
        return 0;

    if (AnimFileInfo->AnimationCount == 0) {
        return 0;
    }
    granny_animation* anim = AnimFileInfo->Animations[0];
    FbxAnimStack* anim_stack = FbxAnimStack::Create(pScene, anim->Name);
    FbxAnimLayer* anim_layer = FbxAnimLayer::Create(pScene, "Layer");

    anim_stack->AddMember(anim_layer);

    for (int i = 0; i < anim->TrackGroupCount; ++i)
        export_animation(pScene, anim, anim->TrackGroups[i], anim_layer);

    return true;
}

bool SaveScene(FbxManager* pManager, FbxDocument* pScene, const char* pFilename, int pFileFormat = -1, bool pEmbedMedia = false)
{
    int lMajor, lMinor, lRevision;
    bool lStatus = true;

    // Create an exporter.
    FbxExporter* lExporter = FbxExporter::Create(pManager, "");

    if (pFileFormat < 0 || pFileFormat >= pManager->GetIOPluginRegistry()->GetWriterFormatCount())
    {
        // Write in fall back format in less no ASCII format found
        pFileFormat = pManager->GetIOPluginRegistry()->GetNativeWriterFormat();

        //Try to export in ASCII if possible
        int lFormatIndex, lFormatCount = pManager->GetIOPluginRegistry()->GetWriterFormatCount();

        for (lFormatIndex = 0; lFormatIndex < lFormatCount; lFormatIndex++)
        {
            if (pManager->GetIOPluginRegistry()->WriterIsFBX(lFormatIndex))
            {
                FbxString lDesc = pManager->GetIOPluginRegistry()->GetWriterFormatDescription(lFormatIndex);
                //const char* lASCII = "ascii";
                const char* lASCII = "binary";
                if (lDesc.Find(lASCII) >= 0)
                {
                    pFileFormat = lFormatIndex;
                    break;
                }
            }
        }
    }

    // Set the export states. By default, the export states are always set to 
    // true except for the option eEXPORT_TEXTURE_AS_EMBEDDED. The code below 
    // shows how to change these states.
    IOS_REF.SetBoolProp(EXP_FBX_MATERIAL, true);
    IOS_REF.SetBoolProp(EXP_FBX_TEXTURE, true);
    IOS_REF.SetBoolProp(EXP_FBX_EMBEDDED, pEmbedMedia);
    IOS_REF.SetBoolProp(EXP_FBX_SHAPE, true);
    IOS_REF.SetBoolProp(EXP_FBX_GOBO, true);
    IOS_REF.SetBoolProp(EXP_FBX_ANIMATION, true);
    IOS_REF.SetBoolProp(EXP_FBX_GLOBAL_SETTINGS, true);

    // Initialize the exporter by providing a filename.
    if (lExporter->Initialize(pFilename, pFileFormat, pManager->GetIOSettings()) == false)
    {
        FBXSDK_printf("Call to FbxExporter::Initialize() failed.\n");
        FBXSDK_printf("Error returned: %s\n\n", lExporter->GetStatus().GetErrorString());
        return false;
    }

    FbxManager::GetFileFormatVersion(lMajor, lMinor, lRevision);
    FBXSDK_printf("FBX file format version %d.%d.%d\n\n", lMajor, lMinor, lRevision);

    // Export the scene.
    lStatus = lExporter->Export(pScene);
    if (!lStatus) {
        FbxString error = lExporter->GetStatus().GetErrorString();
        FBXSDK_printf("Error returned: %s\n\n", error.Buffer());
    }
    else 
    {
        FBXSDK_printf("File %s exported\n", pFilename);
    }
    // Destroy the exporter.
    lExporter->Destroy();
    return lStatus;
}

int main(int argc, char* argv[])
{
    if (argc < 3) {
        cout << "AOE3 gr2 export to FBX\nUsage: gr2fbxexporter.exe <modelfile> <outputfbxfile> [animfile1] [animfile2]\n";
        return 1;
    }
        
    granny_file* GrannyFile = GrannyReadEntireFile(argv[1]);
    if (GrannyFile == 0)
    {
        printf("Error: unable to load %s as a Granny file.\n", argv[1]);
        return EXIT_FAILURE;
    }
    else
    {
        PrintFileStatistics(GrannyFile);
    }

    granny_file_info* FileInfo = GrannyGetFileInfo(GrannyFile);
    if (FileInfo == NULL)
        return 0;


    FbxManager* mSdkManager;
    FbxScene* mScene;
    FbxImporter* mImporter;

    InitializeSdkObjects(mSdkManager, mScene);
    if (mSdkManager)
    {
        CreateScene(mScene, FileInfo);
        AddAnim(mScene, argv[1]);
        for (int i = 3; i < argc; i++) {
            AddAnim(mScene, argv[i]);
        }
        
        bool lResult = SaveScene(mSdkManager, mScene, argv[2]);
        if (lResult == false)
        {
            FBXSDK_printf("\nAn error occurred while saving the scene...\n");
            return 1;
        }
    }

    if (mSdkManager) 
        mSdkManager->Destroy();
    return 0;
}

void PrintFileStatistics(granny_file* GrannyFile)
{
    printf("granny_file\n-----------\n");

    /* Each Granny file contains 1 or more sections that may be loaded, freed, and
       compressed independantly.  By default, the exporter creates several, in order to
       try to segregate ephemeral data that is loaded into the rendering API (for
       instance, vertex buffers and textures) from data that Granny needs to access at
       runtime.  This allows you to release the memory for the ephemeral data if you need
       to save on memory at runtime.  Use the $FreeFileSection API to do this.
     */
    granny_grn_section* SectionArray = GrannyGetGRNSectionArray(GrannyFile->Header);

    printf("File contains: %d sections.\n", GrannyFile->SectionCount);
    {for (granny_int32x Section = 0; Section < GrannyFile->SectionCount; ++Section)
    {
        if (GrannyFile->Sections[Section])
        {
            printf("  Section %d: present", Section);

            // DataSize is the sizeof the data on disk, ExpandedDataSize is the size in memory.  If
            // they match, this section was written without Oodle1 compression.
            if (SectionArray[Section].DataSize ==
                SectionArray[Section].ExpandedDataSize)
            {
                printf(" (uncompressed)\n");
            }
            else
            {
                printf(" (compressed)\n");
            }
        }
        else
        {
            // This will likely never happen, since we haven't called GrannyFreeFileSection
            printf("  Section %d: already freed (or empty)\n", Section);
        }
    }}

    printf("\n");
}

