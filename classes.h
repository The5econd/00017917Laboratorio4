enum lines {NOLINE,SINGLELINE,DOUBLELINE};
enum modes {NOMODE,INT_FLOAT,INT_INT_INT};
enum parameters {ELEMENT_LENGTH,T_RARA,CONST_K,EPSILON,CONST_V,FI,ALFA,SIGMA,MIU};
enum sizes {NODES,ELEMENTS,DIRICHLET};

class item{
protected:
    int id;
    float x;
    int node1;
    int node2;
    float value;
public:
    void setId(int identifier) {
        id = identifier;
    }

    void setX(float x_coord) {
        x = x_coord;
    }

    void setNode1(int node_1) {
        node1 = node_1;
    }

    void setNode2(int node_2) {
        node2 = node_2;
    }

    void setValue(float value_to_assign) {
        value = value_to_assign;
    }

    int getId() {
        return id;
    }

    float getX() {
        return x;
    }

    int getNode1() {
        return node1;
    }

    int getNode2() {
        return node2;
    }

    float getValue() {
        return value;
    }

    virtual void setIntFloat(int n,float r)=0;

    virtual void setIntIntInt(int n1,int n2,int n3)=0;

};

class node: public item{

public:
    void setIntFloat(int identifier, float x_coordinate){
        id = identifier;
        x = x_coordinate;
    }

    void setIntIntInt(int n1,int n2,int n3){
    }

};

class element: public item{

public:
    void setIntFloat(int n1,float r){
    }

    void setIntIntInt(int identifier, int firstnode,int secondnode){
        id = identifier;
        node1 = firstnode;
        node2 = secondnode;
    }

};

class condition: public item{

public:
    void setIntFloat(int node_to_apply, float prescribed_value){
        node1 = node_to_apply;
        value = prescribed_value;
    }

    void setIntIntInt(int n1,int n2,int n3){
    }

};

class mesh{
    float parameters[9];
    int sizes[3];
    node *node_list;
    element *element_list;
    condition *dirichlet_list;
public:
    void setParameters(float l,float t_rara, float k, float epsilon,float v,float fi,float alfa,float sigma,float miu){
        parameters[ELEMENT_LENGTH]=l;
        parameters[T_RARA]=t_rara;
        parameters[CONST_K]=k;
        parameters[EPSILON]=epsilon;
        parameters[CONST_V]=v;
        parameters[FI]=fi;
        parameters[ALFA]=alfa;
        parameters[SIGMA]=sigma;
        parameters[MIU]=miu;
    }
    void setSizes(int nnodes,int neltos,int ndirich){
        sizes[NODES] = nnodes;
        sizes[ELEMENTS] = neltos;
        sizes[DIRICHLET] = ndirich;
    }
    int getSize(int s){
        return sizes[s];
    }
    float getParameter(int p){
        return parameters[p];
    }
    void createData(){
        node_list = new node[sizes[NODES]];
        element_list = new element[sizes[ELEMENTS]];
        dirichlet_list = new condition[sizes[DIRICHLET]];
    }
    node* getNodes(){
        return node_list;
    }
    element* getElements(){
        return element_list;
    }
    condition* getDirichlet(){
        return dirichlet_list;
    }
    node getNode(int i){
        return node_list[i];
    }
    element getElement(int i){
        return element_list[i];
    }
    condition getCondition(int i, int type){
        if(type == DIRICHLET) return dirichlet_list[i];

    }
};
