#include "parse.h"
#include <map>
using namespace std;


namespace {

    // read in dim-dimensional control points into a vector
    vector<Vector3f> readCps(istream &in, unsigned dim)
    {    
        // number of control points    
        unsigned n;
        in >> n;

        cerr << "  " << n << " cps" << endl;
    
        // vector of control points
        vector<Vector3f> cps(n);

        char delim;
        float x;
        float y;
        float z;

        for( unsigned i = 0; i < n; ++i )
        {
            switch (dim)
            {
            case 2:
                in >> delim;
                in >> x;
                in >> y;
                cps[i] = Vector3f( x, y, 0 );
                in >> delim;
                break;
            case 3:
                in >> delim;
                in >> x;
                in >> y;
                in >> z;
                cps[i] = Vector3f( x, y, z );
                in >> delim;
                break;            
            default:
                abort();
            }
        }

        return cps;
    }
}

Curve evalCatMullRom(const vector< Vector3f >& P, unsigned steps)
{
    // Check
    if (P.size() < 4)
    {
        cerr << "evalCatMullRom must be called with 4 or more control points." << endl;
        exit(0);
    }

    // TODO:
    // Implement CatMull Rom splines that connect lines through the control points

    cerr << "\t>>> evalCatMullRom has been called with the following input:" << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): " << endl;
    for (unsigned i = 0; i < P.size(); ++i)
    {
        cerr << "\t>>> " << P[i] << endl;
    }

    cerr << "\t>>> Steps (type steps): " << steps << endl;
    //cerr << "\t>>> Returning empty curve." << endl;

    int seg;			//Determine number of segments if 4 or 
    //more control points
    if (P.size() == 4) {
        seg = 1;
    }
    else {
        seg = (P.size() - 3);
    }

    Curve R((seg) * (steps)+1);

    Vector3f Binit;		//Initial Binormal (arbitrary only at beginning)
    Vector3f Blast;		//Most recent Binormal

    bool beginning = true;

    for (unsigned i = 0; i < P.size() - 3; ++i)
    {
        if (beginning) {
            Binit = Vector3f(0.f, 0.f, 1.f);
        }
        else {
            Binit = Blast;
        }

        for (unsigned delta = 0; delta <= steps; ++delta)
        {
            float t = float(delta) / steps;

            //Point Matrix
            Matrix4f relativePoints(P[i + 0][0], P[i + 1][0], P[i + 2][0], P[i + 3][0],
                P[i + 0][1], P[i + 1][1], P[i + 2][1], P[i + 3][1],
                P[i + 0][2], P[i + 1][2], P[i + 2][2], P[i + 3][2],
                0.f, 0.f, 0.f, 0.f);
            //V Matrix
            Matrix4f MatV(0.f / 2, -1.f / 2, 2.f / 2, -1.f / 2,
                2.f / 2, 0.f, -5.f / 2, 3.f / 2,
                0.f / 2, 1.f / 2, 4.f / 2, -3.f / 2,
                0.f, 0.f, -1.f / 2, 1.f / 2);
            //T Matrix
            Matrix4f MatT(-1.f / 2, 4.f / 2, -3.f / 2, 0.f / 2,
                0.f / 2, -10.f / 2, 9.f / 2, 0.f / 2,
                1.f / 2, 8.f / 2, -9.f / 2, 0.f / 2,
                0.f, -2.f / 2, 3.f / 2, 0.f / 2);

            //polynomial (t) Matrix
            Vector4f timeT(1, t, t * t, t * t * t);

            //Calculate V Vector
            R[(i * steps) + delta].V = Vector3f((relativePoints * MatV * timeT)[0],
                (relativePoints * MatV * timeT)[1],
                (relativePoints * MatV * timeT)[2]);
            //Calculate Tangent
            R[(i * steps) + delta].T = Vector3f((relativePoints * MatT * timeT)[0],
                (relativePoints * MatT * timeT)[1],
                (relativePoints * MatT * timeT)[2]).normalized();
            //Calculate Normal
            R[(i * steps) + delta].N = Vector3f::cross(Binit,
                R[(i * steps) + delta].T).normalized();
            //Calulate Binormal
            R[(i * steps) + delta].B = Vector3f::cross(R[(i * steps) + delta].T,
                R[(i * steps) + delta].N).normalized();
            //Keep track of current Binormal
            Binit = R[(i * steps) + delta].B;

            beginning = false;
            Blast = Binit;
        }
    }

    // Return an empty curve right now.
    //return Curve();
    return R;
}

bool parseFile(istream &in,
               vector<vector<Vector3f> > &ctrlPoints, 
               vector<Curve>             &curves,
               vector<string>            &curveNames,
               vector<Surface>           &surfaces,
               vector<string>            &surfaceNames)
{
    ctrlPoints.clear();
    curves.clear();
    curveNames.clear();
    surfaces.clear();
    surfaceNames.clear();    
    
    string objType;

    // For looking up curve indices by name
    map<string,unsigned> curveIndex;

    // For looking up surface indices by name
    map<string,unsigned> surfaceIndex;
        
    // For storing dimension of curve
    vector<unsigned> dims;

    unsigned counter = 0;
    
    while (in >> objType) 
    {
        cerr << ">object " << counter++ << endl;
        string objName;
        in >> objName;

        bool named = (objName != ".");
        
        vector<Vector3f> cpsToAdd;
        
        if (curveIndex.find(objName) != curveIndex.end() ||
            surfaceIndex.find(objName) != surfaceIndex.end())
        {
            cerr << "error, [" << objName << "] already exists" << endl;
            return false;
        }

        unsigned steps;

        if (objType == "bez2")
        {
            in >> steps;
            cerr << " reading bez2 " << "[" << objName << "]" << endl;
            curves.push_back( evalBezier(cpsToAdd = readCps(in, 2), steps) );
            curveNames.push_back(objName);
            dims.push_back(2);
            if (named) curveIndex[objName] = dims.size()-1;
            
        }
        else if (objType == "bsp2")
        {
            cerr << " reading bsp2 " << "[" << objName << "]" << endl;
            in >> steps;
            curves.push_back( evalBspline(cpsToAdd = readCps(in, 2), steps) );
            curveNames.push_back(objName);
            dims.push_back(2);
            if (named) curveIndex[objName] = dims.size()-1;
        }

        else if (objType == "bez3")
        {
            cerr << " reading bez3 " << "[" << objName << "]" << endl;
            in >> steps;
            curves.push_back( evalBezier(cpsToAdd = readCps(in, 3), steps) );
            curveNames.push_back(objName);
            dims.push_back(3);
            if (named) curveIndex[objName] = dims.size()-1;

        }
        else if (objType == "bsp3")
        {
            cerr << " reading bsp3 " << "[" << objName << "]" << endl;
            in >> steps;
            curves.push_back( evalBspline(cpsToAdd = readCps(in, 3), steps) );
            curveNames.push_back(objName);
            dims.push_back(3);
            if (named) curveIndex[objName] = dims.size()-1;
        }


//CatMull Rom splines


        else if (objType == "cmr2")
        {
            cerr << " reading cmr2 " << "[" << objName << "]" << endl;
            in >> steps;
            curves.push_back(evalCatMullRom(cpsToAdd = readCps(in, 2), steps) );
            curveNames.push_back(objName);
            dims.push_back(2);
            if (named) curveIndex[objName] = dims.size()-1;
        }
	else if (objType == "cmr3")
        {
            cerr << " reading cmr3 " << "[" << objName << "]" << endl;
            in >> steps;
            curves.push_back(evalCatMullRom(cpsToAdd = readCps(in, 3), steps) );
            curveNames.push_back(objName);
            dims.push_back(3);
            if (named) curveIndex[objName] = dims.size()-1;
        }




        else if (objType == "srev")
        {
            cerr << " reading srev " << "[" << objName << "]" << endl;
            in >> steps;

            // Name of the profile curve
            string profName;
            in >> profName;

            cerr << "  profile [" << profName << "]" << endl;
            
            map<string,unsigned>::const_iterator it = curveIndex.find(profName);

            // Failure checks
            if (it == curveIndex.end()) {                
                cerr << "failed: [" << profName << "] doesn't exist!" << endl; return false;
            }
            if (dims[it->second] != 2) {
                cerr << "failed: [" << profName << "] isn't 2d!" << endl; return false;
            }

            // Make the surface
            surfaces.push_back( makeSurfRev( curves[it->second], steps ) );
            surfaceNames.push_back(objName);
            if (named) surfaceIndex[objName] = surfaceNames.size()-1;
        }
        else if (objType == "gcyl")
        {
            cerr << " reading gcyl " << "[" << objName << "]" << endl;
            
            // Name of the profile curve and sweep curve
            string profName, sweepName;
            in >> profName >> sweepName;

            cerr << "  profile [" << profName << "], sweep [" << sweepName << "]" << endl;

            map<string,unsigned>::const_iterator itP, itS;

            // Failure checks for profile
            itP = curveIndex.find(profName);
            
            if (itP == curveIndex.end()) {                
                cerr << "failed: [" << profName << "] doesn't exist!" << endl; return false;
            }
            if (dims[itP->second] != 2) {
                cerr << "failed: [" << profName << "] isn't 2d!" << endl; return false;
            }

            // Failure checks for sweep
            itS = curveIndex.find(sweepName);
            if (itS == curveIndex.end()) {                
                cerr << "failed: [" << sweepName << "] doesn't exist!" << endl; return false;
            }

            // Make the surface
            surfaces.push_back( makeGenCyl( curves[itP->second], curves[itS->second] ) );
            surfaceNames.push_back(objName);
            if (named) surfaceIndex[objName] = surfaceNames.size()-1;

        }
        else if (objType == "circ")
        {
            cerr << " reading circ " << "[" << objName << "]" << endl;

            unsigned steps;
            float rad;
            in >> steps >> rad;
            cerr << "  radius [" << rad << "]" << endl;

            curves.push_back( evalCircle(rad, steps) );
            curveNames.push_back(objName);
            dims.push_back(2);
            if (named) curveIndex[objName] = dims.size()-1;
        }
        else
        {
            cerr << "failed: type " << objType << " unrecognized." << endl;
            return false;
        }

        ctrlPoints.push_back(cpsToAdd);
    }

    return true;
}




