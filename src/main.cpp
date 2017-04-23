#include <iostream>
#include "SimpleGraph.h"
#include <fstream>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <cmath>
#include <ctime>
using namespace std;

/**
 * cleans up main
 * prints given welcome message for user
 */
void Welcome() {
    cout << "Welcome to CS106L GraphViz!" << endl;
    cout << "This program uses a force-directed graph layout algorithm" << endl;
    cout << "to render sleek, snazzy pictures of various graphs." << endl;
    cout << endl;
}

/**
 * I keep track of dx's and dy's using a vector
 * of these structs
 */
struct Force {
    //even though the SimpleGraph's Node object also consists of  two doubles,
    //I wanted to name a struct separately to not get confused
    double dx;
    double dy;
};

/**
 * Takes a line of input, and writes the two edges by references
 * No error checking - good input is guaranteed
 * @brief parse_ints
 * @param line
 * @param edge1
 * @param edge2
 */
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

/**
 * For testing purposes - prints single node
 * @brief print_node
 * @param n
 */
void print_node(const Node &n) {
    cout << "(" << n.x << ", " << n.y << ")" << endl;
}

/**
 * For testing purposes - prints vector of node
 * @brief print_nodes
 * @param nvec
 */
void print_nodes(const vector<struct Node> &nvec) {
    for (unsigned int i=0; i<nvec.size(); i++) {
        print_node(nvec[i]);
    }
    cout << endl;
}

/**
 * For testing - prints a single force struct
 * @brief print_force
 * @param f
 */
void print_force(const Force &f) {
    cout << "[" << f.dx << ", " << f.dy << "]" << endl;
}

/**
 * For testing - prints vector of forces
 * @brief print_forces
 * @param force_vec
 */

void print_forces(const vector<struct Force> &force_vec) {
    for (unsigned int i=0; i<force_vec.size(); i++) {
        print_force(force_vec[i]);
    }
    cout << endl;
}

/**
 * For testing - prints a single edge
 * @brief print_edge
 * @param e
 */
void print_edge(const Edge &e) {
    cout << e.start << " -> " << e.end << endl;
}

/**
 * For testing - prints vector of edges
 * @brief print_edges
 * @param evec
 */
void print_edges(const vector<Edge> &evec) {
    for (unsigned int i=0; i<evec.size(); i++) {
        print_edge(evec[i]);
    }
}

/**
 * Returns theta between two nodes
 * Calculation is same for attraction and repulsion
 * @brief calc_theta
 * @param n0
 * @param n1
 * @return
 */
double calc_theta(const Node &n0, const Node &n1) {
    return atan2(n1.y - n0.y, n1.x - n0.x);
}

/**
  Returns constant F-repel
 * @brief f_rep
 * @param n0
 * @param n1
 * @param k
 * @return
 */
double f_rep(const Node &n0, const Node &n1, const double &k) {
    double ysqr = pow((n1.y - n0.y),2.0);
    ysqr += pow((n1.x - n0.x), 2.0);
    return (k / sqrt(ysqr));
}



/**
 * Returns constant F-attraction
 * @brief f_att
 * @param n0
 * @param n1
 * @param k
 * @return
 */
double f_att(const Node &n0, const Node &n1, const double &k) {
    double ysqr = pow((n1.y - n0.y),2.0);
    ysqr += pow((n1.x - n0.x), 2.0);
    return (k * ysqr);
}


/**
 * Adjusts for attraction between every edge
 * Stores results in vector of Forces
 * @brief do_update1
 * @param sg
 * @param forces
 */
void do_update1(const SimpleGraph &sg, vector<struct Force> &forces) {
    double ka = 0.001;
    Node n0, n1;
    Edge e;
    for (unsigned int i=0; i<sg.edges.size(); i++) {
        e = sg.edges[i];
        n0 = sg.nodes[e.start];
        n1 = sg.nodes[e.end];
        double theta = calc_theta(n0, n1);
        double fatt = f_att(n0, n1, ka);
        double xval2 = fatt * cos(theta);
        double yval2 = fatt * sin(theta);

        forces[e.start].dx += xval2;
        forces[e.start].dy += yval2;
        forces[e.end].dx -= xval2;
        forces[e.end].dy -= yval2;
    }
}

/**
 * Adjusts for repulsion between all pairs of nodes
 * Stores results in a vector of Forces
 * @brief do_update2
 * @param nvec
 * @param forces
 */
void do_update2(const vector<struct Node> &nvec, vector<struct Force> &forces) {
    double kr = 0.001;
    Node n0, n1;
    double f, theta;
    double xval, yval;
    for (unsigned int i=0; i<nvec.size(); i++) {
        for (unsigned int j=i+1; j<nvec.size(); j++) {
            //cout << i << " , " << j << endl;
            n0 = nvec[i];
            n1 = nvec[j];
            f = f_rep(n0, n1, kr);
            theta = calc_theta(n0, n1);
            xval = f * cos(theta);
            yval = f * sin(theta);

            forces[i].dx -= xval;
            forces[i].dy -= yval;
            forces[j].dx += xval;
            forces[j].dy += yval;
        }
    }
}

/**
 * Changes coordinates of SG's nodes based on Forces
 * @brief make_change
 * @param sg
 * @param forces
 */
void make_change(SimpleGraph &sg, vector<struct Force> &forces) {

    for (unsigned int i=0; i<forces.size(); i++) {
        Force f = forces[i];
        sg.nodes[i].x += f.dx;
        sg.nodes[i].y += f.dy;
    }
}

/**
 * Creates vector of Forces - one force for each node
 * Dx and Dy for each node = 0
 * @brief init_forces
 * @param forces
 * @param num_nodes
 */
void init_forces(vector<struct Force> &forces, const int num_nodes) {
    struct Force null_force;
    null_force.dx = 0;
    null_force.dy = 0;

    for (int i=0; i<num_nodes; i++) {
        forces.push_back(null_force);
    }
}

/**
 * Sets each Force's Dx and Dy = 0
 * @brief reset_forces
 * @param forces
 * @param num_nodes
 */
void reset_forces(vector<struct Force> &forces, const int num_nodes) {
    for (int i=0; i<num_nodes; i++) {
        forces[i].dx = 0;
        forces[i].dy = 0;
    }
}

/**
 * Gets filename from user; reads file into graph
 * @brief read_file
 * @param Pi
 * @param sg
 */
void read_file(const double &Pi, SimpleGraph &sg) {

    bool valid = false;

    string filename;
    int fd;
    cout << "To begin, enter the name of your graph file." << endl;
    int numtries = 0;
    //cout << "num tries = " << numtries << endl;
    while (valid == false) {
      //  cout << "num tries = " << numtries << endl;
        if (numtries > 1) cout << "That file was not found. Try again: " << endl;
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

    for (int i=num_nodes; i>0; i--) {
        Node n;
        double xc = cos(((2*Pi*i) / num_nodes));
        double yc = sin(((2*Pi*i) / num_nodes));
        n.x = xc;
        n.y = yc;
        sg.nodes.push_back(n);
    }

    int edge1, edge2;
    int num_edges = 0;
    for (int i=svec.size()-1; i>0; i--) {
        parse_ints(svec[i], edge1, edge2);
        Edge e;
        e.start = edge1;
        e.end = edge2;
        sg.edges.push_back(e);
        num_edges ++;
    }
}

/**
 * Did not end up using this function (yet)
 * @brief isnum
 * @param str
 * @return
 */
bool isnum(string &str) {
    for (unsigned int i=0; i<str.size(); i++) {
        if (!isdigit(str[i]) && str[i] != '.') {
            return false;
        }
    }
    return true;
}

/**
 * Gets time from user.
 * This problem needs additional error checking
 * If the user enters a string (non-numeric) of length
 * two or greater, then the program basically fails.
 * I tried for a long time to fix this (using readline, using
 * stringstreams, but was confounded at every turn).
 * Nevertheless, if a normal double is entered, it works fine.
 * @brief get_time
 * @param result
 */
void get_time(double &result) {
    bool valid = false;
    int numtries = 0;
    long double t;
    while (valid == false) {
        if (numtries > 0) cout << "That's not a valid length of time. Try again:" << endl;
        cout << "How many seconds: "; cin >> t;
        if (t > 0) valid = true;
        numtries ++;
        cin.clear();
    }
    result = t;
}

int main() {
    Welcome();
    SimpleGraph sg;

    const double Pi = 3.14159265358979323;
    double d;
    get_time(d);
//    cout << "will run for: " << d << " seconds." << endl;
    read_file(Pi, sg);
    InitGraphVisualizer(sg);
    int num_nodes = sg.nodes.size();
    DrawGraph(sg);
    vector<struct Force> forces;
    init_forces(forces, num_nodes);
    vector<struct Node> nodez;

    time_t startTime = time(NULL);
    double elapsedTime = 0;
    while(elapsedTime < d) {
        nodez = sg.nodes;
        do_update1(sg, forces);
        do_update2(nodez, forces);
        make_change(sg, forces);
        DrawGraph(sg);
        reset_forces(forces, num_nodes);
        elapsedTime = difftime(time(NULL), startTime);
    }
    cout << "Thanks for playing!" << endl;
    return 0;
}
