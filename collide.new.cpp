#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <tipsy_wrapper.hh>
#include <particle.hh>

using std::cout;
using std::endl;
using std::cin;

class collide
{
private:
    double x1;
    double y1;
    double z1;

    double x2;
    double y2;
    double z2;

    double vx1;
    double vy1;
    double vz1;
    double vx2;
    double vy2;
    double vz2;

public:
    std::string profile1;
    std::string profile2;
    std::string output;

    std::vector<Particle> p;

    collide(std::string in1, std::string in2, std::string out, double dx1, double dy1, double dz1,
            double dx2, double dy2, double dz2, double v_x1, double v_y1, double v_z1, double v_x2,
            double v_y2, double v_z2)
    {
        // Initialise variables
        profile1 = in1;
        profile2 = in2;
        output = out;

        x1 = dx1;
        y1 = dy1;
        z1 = dz1;

        x2 = dx2;
        y2 = dy2;
        z2 = dz2;

        vx1 = v_x1;
        vy1 = v_y1;
        vz1 = v_z1;

        vx2 = v_x2;
        vy2 = v_y2;
        vz2 = v_z2;

        // Later we want to do calculations in the COM system but for now we assume that M1 = M2
        // and displace both bodies in a distance of 10 R_E along the x-axis

        // At the moment the user can only choose the displacement along the y-axis
        // and the speed in x direction

    }

    double get_R(std::string file)
    {
        double R, r;
        TipsyInFile in;

        R = 0.0;

        in.open(file);

        std::vector<Particle> p(in.size());

        for (size_t i = 0; i < p.size(); ++i)
        {
            // Read one particle from Tipsy file
            in >> p[i];

            r = sqrt(p[i].x*p[i].x + p[i].y*p[i].y + p[i].z*p[i].z);
            if (r > R)
            {
                R = r;
            }
        }

        return R;
    }

    void write_file()
    {
        TipsyInFile in;
        Particle particle;

        // Read Profile1
        cout << "Reading input file " << profile1 << endl;
        in.open(profile1);

        cout << "N: " << in.size() << " t: " << in.time() << endl;

        while (!in.eof())
        {
            in >> particle;

            // Change positions
            particle.x = particle.x + x1;
            particle.y = particle.y + y1;
            particle.z = particle.z + z1;

            // Change velocities
            particle.vx += vx1;
            particle.vy += vy1;
            particle.vz += vz1;

            p.push_back(particle);
        }
        in.close();

        // Read Profile2
        cout << "Reading input file " << profile2 << endl;
        in.open(profile2);

        cout << "N: " << in.size() << " t: " << in.time() << endl;

        while (!in.eof())
        {
            in >> particle;

            // Change positions
            particle.x = particle.x + x2;
            particle.y = particle.y + y2;
            particle.z = particle.z + z2;

            // Change velocities
            particle.vx += vx2;
            particle.vy += vy2;
            particle.vz += vz2;

            p.push_back(particle);
        }
        in.close();

        // Write output file
        TipsyOutFile out;
        out.open(output);
        out.time(0.0);

        for (size_t i = 0; i < p.size(); ++i)
        {
            out << p[i];
        }

        cout << "Writing output file " << output << endl;
        cout << "N: " << p.size();

        out.close();

        cout << "Done." << endl;

    }
};

double get_R(std::string file)
{
    double R, r;
    TipsyInFile in;

    R = 0.0;

    in.open(file);

    std::vector<Particle> p(in.size());

    for (size_t i = 0; i < p.size(); ++i)
    {
        // Read one particle from Tipsy file
        in >> p[i];

        r = sqrt(p[i].x*p[i].x + p[i].y*p[i].y + p[i].z*p[i].z);
        if (r > R)
        {
            R = r;
        }
    }

    return R;
}
double get_M(std::string file)
{
    double M;
    TipsyInFile in;

    M = 0.0;

    in.open(file);

    std::vector<Particle> p(in.size());

    for (size_t i = 0; i < p.size(); ++i)
    {
        // Read one particle from Tipsy file
        in >> p[i];

        M = M + p[i].m;
    }

    return M;
}
int main(int argc, char** argv)
{    
    double R1, R2, M1, M2, b, gamma;
    double x1, y1, x2, y2, deltax, deltay, v1, v2, deltav;

    std::string in1;
    std::string in2;
    std::string out;

    R1 = 0.0;
    R2 = 0.0;
    M1 = 0.0;
    M2 = 0.0;

    b = 0.0;
    gamma = 0.0;

    x1 = 0;
    y1 = 0;
    x2 = 0;
    y2 = 0;
    deltax = 0;
    deltay = 0;

    if (argc != 8) {
        cout << "Usage: collide <b> <v> <Rt> <Ri> <profile1> <profile2> <output.std>" << endl;
        return 1;
    }

    b = atof(argv[1]);
    assert(b >= 0 && "Impact parameter has to be larger than zero");
    deltav = atof(argv[2]);
    assert(deltav >= 0 && "Impact velocity has to be larger than zero");
	
	// Read radius from the command line
    R1 = atof(argv[3]);
    assert(R1 >= 0 && "R1 has to be larger than zero");
    R2 = atof(argv[4]);
    assert(R2 >= 0 && "R2 has to be larger than zero");
    in1 = argv[5];
    in2 = argv[6];
    out = argv[7];

    /* Determine total mass and radius of the equilibrium models */

//    R1 = get_R(in1);
//    R2 = get_R(in2);

    M1 = get_M(in1);
    M2 = get_M(in2);

    /* Determine the coordinates in the COM system */

    gamma = M2/(M1 + M2);

//    deltax = 5.0;
    deltax = (R1 + R2)*1.1;

    deltay = (R1 + R2)*b;

    x1 = -deltax*gamma;
    x2 = deltax*(1 - gamma);
    y1 = -deltay*gamma;
    y2 = deltay*(1 - gamma);

    v1 = deltav*gamma;
    v2 = -deltav*(1 - gamma);

    cout << "Profile1:" << in2 << " Profile2: " << in2 << "Ouput: " << out << " b: " << b << " v: " << deltav << endl;

    cout << "M1: " << M1 << " M2: " << M2 << " gamma: " << gamma << endl;
    cout << "x1: " << x1 << " y1: " << y1 << " x2: " << x2 << " y2: " << y2 << endl;
    cout << "v1: " << v1 << " v2: " << v2 << " R1: " << R1 << " R2: " << R2 << endl;

    collide c(in1, in2, out, x1, y1, 0.0, x2, y2, 0.0, v1, 0.0, 0.0, v2, 0.0, 0.0);
/*
        try {
            p.read_datafile();
        }
        catch (const std::ifstream::failure& e)
        {
            cout << "Error while opening/reading file: " << e.what() << endl;
            return 1;
        }
*/
        c.write_file();

        return 0;
}
