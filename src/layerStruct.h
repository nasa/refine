typedef struct Normal Normal;
struct Normal {
  int constrained;
  int root, tip;
  double direction[3];
  double height;
  double initialheight;
  double length;
  double maxlength;
  double rate;
  bool terminated;
};

typedef struct Triangle Triangle;
struct Triangle {
  int normal[3];
  int constrainedSide[3];
  int parentGeomEdge[3];
};

typedef struct Blend Blend;
struct Blend {
  int nodes[2];
  int normal[4];
  int edgeId[2];
  int oldnormal[4];
};

struct Layer {
  Grid *grid;
  int maxtriangle, ntriangle;
  Triangle *triangle;
  int maxParentGeomFace, nParentGeomFace, *ParentGeomFace;
  int maxblend, nblend;
  Blend *blend;
  int maxnormal, nnormal;
  Normal *normal;
  int *globalNode2Normal;
  int nConstrainingGeometry, *constrainingGeometry;
  Adj *adj;
  Near *nearTree;
  bool mixedElementMode;

  bool *cellInLayer;
  bool *faceInLayer;
  bool *edgeInLayer;

  FILE *tecplotFile;
};
