void do_repulse(vector<struct Force> &forces, const vector<struct Edge> &edges, const vector<struct Node> &nodes) {
    double KRepel = 0.001;
   // print_edges(edges);
   // print_nodes(nodes);
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

void do_attract(vector<struct Force> &forces, const vector<struct Edge> &edges, const vector<struct Node> &nodes) {
    double KAttract = 0.001;
    //print_edges(edges);
    //print_nodes(nodes);
    //print_forces(forces);
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
    //update_graph(forces, sg);
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
