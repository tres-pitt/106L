#include <iostream>
#include "SimpleGraph.h"
#include <fstream>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <cmath>

using namespace std;

void Welcome() {
    cout << "Welcome to CS106L GraphViz!" << endl;
    cout << "This program uses a force-directed graph layout algorithm" << endl;
    cout << "to render sleek, snazzy pictures of various graphs." << endl;
    cout << endl;
}

//even though the SimpleGraph's Node object also has two doubles,
//I wanted to name a struct separately to not get confused
struct Force {
    double dx;
    double dy;
};

void parse_ints(string &line, int &edge1, int &edge2) {
    char ch = line[0];
    if (isdigit(ch)==false) {
        cout << "parse ints error" << endl;
        return;
    } else {
        int counter = 0;
        bool found = false;
        while (found == false) {
            counter ++;
            char c = line[counter];
            if (isdigit(c)==false) {
                found=true;
            }
        }
        edge1 = stoi(line.substr(0, counter));
        edge2 = stoi(line.substr(counter, line.size()-counter));
    }
}


void print_node(const Node n) {
    cout << "(" << n.x << ", " << n.y << ")" << endl;
}

void print_nodes(const vector<struct Node> &nvec) {
    for (unsigned int i=0; i<nvec.size(); i++) {
        print_node(nvec[i]);
    }
    cout << endl;
}

void print_force(const Force f) {
    cout << "[" << f.dx << ", " << f.dy << "]" << endl;
}

void print_forces(const vector<struct Force> &force_vec) {
    for (unsigned int i=0; i<force_vec.size(); i++) {
        print_force(force_vec[i]);
    }
    cout << endl;
}

void print_edge(const Edge &e) {
    cout << e.start << " -> " << e.end << endl;
}

void print_edges(const vector<Edge> evec) {
    for (unsigned int i=0; i<evec.size(); i++) {
        print_edge(evec[i]);
    }
}

double calc_theta(const Node &n0, const Node &n1) {
    return atan2(n1.y - n0.y, n1.x - n0.x);
}

double f_rep(const Node &n0, const Node &n1, const double &k) {
    double ysqr = pow((n1.y - n0.y),2.0);
    ysqr += pow((n1.x - n0.x), 2.0);
    return (k / sqrt(ysqr));
}


void do_repulse(vector<struct Force> &forces, const vector<struct Edge> &edges, const vector<struct Node> &nodes) {
    double KRepel = 0.001;
    print_edges(edges);
    print_nodes(nodes);
    for (unsigned int i=0; i<edges.size(); i++) {
        struct Edge e = edges[i];
        Node n0 = nodes[e.start];
        Node n1 = nodes[e.end];
        double theta = calc_theta(n0, n1);
        double f = f_rep(n0, n1, KRepel);
        double xval = f * cos(theta);
        double yval = f * sin(theta);

        Force f0 = forces[e.start];
        Force f1 = forces[e.end];
        double f0x = f0.dx - xval;
        double f0y = f0.dy - yval;
        double f1x = f1.dx + xval;
        double f1y = f1.dy + yval;

        f0.dx = f0x;
        f0.dy = f0y;
        f1.dx = f1x;
        f1.dy = f1y;

        forces[e.start] = f0;
        forces[e.end] = f1;
    }
}

double f_att(const Node &n0, const Node &n1, const double &k) {
    double ysqr = pow((n1.y - n0.y),2.0);
    ysqr += pow((n1.x - n0.x), 2.0);
    return (k * ysqr);
}

void do_attract(vector<struct Force> &forces, const vector<struct Edge> &edges, const vector<struct Node> &nodes) {
    double KAttract = 0.001;
    print_edges(edges);
    print_nodes(nodes);
    print_forces(forces);
    for (unsigned int i=0; i<edges.size(); i++) {
        struct Edge e = edges[i];
        Node n0 = nodes[e.start];
        Node n1 = nodes[e.end];
        double theta = calc_theta(n0, n1);
        double f = f_att(n0, n1, KAttract);
        double xval = f * cos(theta);
        double yval = f * sin(theta);

        Force f0 = forces[e.start];
        Force f1 = forces[e.end];
        double f0x = f0.dx + xval;
        double f0y = f0.dy + yval;
        double f1x = f1.dx - xval;
        double f1y = f1.dy - yval;

        f0.dx = f0x;
        f0.dy = f0y;
        f1.dx = f1x;
        f1.dy = f1y;

        forces[e.start] = f0;
        forces[e.end] = f1;
    }
}

void update_graph(vector<struct Force> &forces, SimpleGraph &sg) {
    vector<struct Node> nodez = sg.nodes;

    for (unsigned int i=0; i<nodez.size(); i++) {
        sg.nodes[i].x += forces[i].dx;
        sg.nodes[i].y += forces[i].dy;
//        Node n = nodez[i];
        //Force f = forces[i];
//        double newx = n.x + f.dx;
//        double newy = n.y + f.dy;
//        n.x = newx;
//        n.y = newy;
//        nodez[i] = n;
    }
    //sg.nodes = nodez;
}

void init_forces(vector<struct Force> &forces, const int num_nodes) {
    struct Force null_force;
    null_force.dx = 0;
    null_force.dy = 0;
    for (int i=0; i<num_nodes; i++) {
        forces.push_back(null_force);
    }
}

void calc_forces(vector<struct Force> &forces, SimpleGraph &sg, const int num_nodes) {
    init_forces(forces, num_nodes);
    vector<struct Edge> edges = sg.edges;
    vector<struct Node> nodes = sg.nodes;
//    print_forces(forces);
//    cout << endl;
    do_repulse(forces, edges, nodes);
//    print_forces(forces);
//    cout << endl;
    do_attract(forces, edges, nodes);
//    print_forces(forces);
//    cout << endl;
//    print_nodes(sg.nodes);
    update_graph(forces, sg);
//    print_nodes(sg.nodes);
    DrawGraph(sg);
}

void calculate_forces(vector<struct Force> &forces, SimpleGraph &sg) {
    double KRepel = 0.001;
    double KAttract = 0.001;
    vector<Node> curr_nodes = sg.nodes;
    vector<Edge> curr_edges = sg.edges;

    for (unsigned int i=0; i<curr_edges.size(); i++) {

        Edge e = curr_edges[i];
        int ind1 = e.start;
        int ind2 = e.end;
        cout << "ind1 = " << ind1 << ", ind2 = " << ind2 << endl;
        Node n0 = curr_nodes[ind1];
        Node n1 = curr_nodes[ind2];
        print_node(n0);
        print_node(n1);
        //vector<struct Force> fvec1 = forces[e.start];
        //vector<struct Force> fvec2 = forces[e.end];

        //may need to call new operator here - ?
        Force f0 = forces[ind1];
        Force f1 = forces[ind2];
        cout << "Krepel = " << KRepel << endl;
        double KRepel = 0.001;
        double dist = pow((n1.y - n0.y), 2)  + pow((n1.x - n0.x), 2);
        double frepel = KRepel / sqrt(dist);
        double theta = atan2(n1.y - n0.y, n1.x - n0.x);

        cout << "dist = " << dist << endl;
        cout << "frepel = " << frepel << endl;
        cout << "theta = " << theta << endl;

        double dx = frepel * cos(theta);
        double dy = frepel * sin(theta);
        double dx_new0 = f0.dx - dx;
        double dy_new0 = f0.dy - dy;
        double dx_new1 = f1.dx + dx;
        double dy_new1 = f1.dy + dy;

        double fattract = KAttract * dist;
        double dx2 = fattract * cos(theta);
        double dy2 = fattract * sin(theta);

        dx_new0 += dx2;
        dy_new0 += dy2;
        dx_new1 -= dx2;
        dy_new1 -= dy2;

        f0.dx = dx_new0;
        f0.dy = dy_new0;
        f1.dx = dx_new1;
        f1.dy = dy_new1;
        forces[ind1] = f0;
        forces[ind2] = f1;
        print_force(f0);
        print_force(f1);
        cout << endl;
    }
}


void fix_pos(vector<struct Force> &forces, SimpleGraph &sg) {
    cout << "begin fix_pos" << endl;
    vector<Node> nodes = sg.nodes;
    vector<Node> newnodes;
    for (unsigned int i=0; i<nodes.size(); i++) {
        Node n = nodes[i];
        print_node(n);
        double curr_x = n.x;
        double curr_y = n.y;

        Force f = forces[i];
        print_force(f);

        curr_x += f.dx;
        curr_y += f.dy;
        Node nn;
        nn.x = curr_x;
        nn.y = curr_y;
        cout << "{" << nn.x << ", " << nn.y << "}" << endl;
        //newnodes[i] = nn;
        newnodes.push_back(nn);
        cout << "end it" << endl;
    }
    sg.nodes = newnodes;
    cout << "newnodes: " << endl;
    for (unsigned int i=0; i<newnodes.size(); i++) {
        Node n = newnodes[i];
        print_node(n);
    }
}


int main() {
    const double Pi = 3.14159265358979323;

    Welcome();
    bool valid = false;
    int numtries = 0;
    string filename;
    int fd;
    cout << "To begin, enter the name of your graph file." << endl;
    while (valid == false) {
        if (numtries > 0) cout << "That file was not found. Try again." << endl;
        getline(cin, filename);
        fd = open(filename.c_str(), O_RDONLY);
        if (fd != -1) valid = true;
        numtries ++;
    }
    string line;
    ifstream fs(filename.c_str());
    vector<string> svec;
    if (fs.is_open()) {
        while (getline(fs, line)) {
            svec.push_back(line);
        }
    } else cout << "unable to open file" << endl;
    int num_nodes = stoi((svec[0]));
    cout << "num_nodes = " << num_nodes << endl;
    SimpleGraph sg;
    InitGraphVisualizer(sg);
  //  vector<struct Node> nvec;
//    for (int i=num_nodes; i>0; i--) {
//    for (int i=0; i<num_nodes; i++) {
//        sg.nodes[i].x = cos(((2*Pi*i) / num_nodes));
//        sg.nodes[i].y = sin(((2*Pi*i) / num_nodes));
//    }

    for (int i=0; i<num_nodes; i++) {
        Node n;
        double xc = cos(((2*Pi*i) / num_nodes));
        double yc = sin(((2*Pi*i) / num_nodes));
        n.x = xc;
        n.y = yc;
        sg.nodes.push_back(n);
        //nvec[i] = n;
    }
    //print_nodes(sg.nodes);
    //vector<struct Edge> evec;
    cout << "here" << endl;
    int edge1, edge2;
    int num_edges = 0;
//    for (int i=svec.size()-1; i>0; i--) {
    for (unsigned int i=0; i<svec.size()-1; i++) {
        parse_ints(svec[i], edge1, edge2);
        Edge e;
        e.start = edge1;
        e.end = edge2;
        sg.edges.push_back(e);
     //   evec[i] = e;
        num_edges ++;
    }
    //print_edges(sg.edges);
    //sg.nodes = nvec;
    //sg.edges = evec;
    DrawGraph(sg);
//    vector<struct Force> forces;
//    while(true) {
//        calc_forces(forces, sg, num_nodes);
//    }
    return 0;
}
