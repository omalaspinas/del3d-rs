int init_hull3D_compact(std::vector<R3> &pts, std::vector<Tri> &hull)
{
    int nump = (int)pts.size();
    std::vector<Snork> norts;
    hull.reserve(nump * 2);

    // keep track of covered (dead) triangles.
    std::vector<int> Dlist;
    for (int d = 0; d < 64; d++)
    {
        Dlist.push_back(-1);
    }
    int numD = 64;  // number of dead triangle slots.
    int D_idx = -1; // index to next dead triangle slot to be retuned to use.

    // keep track of last triangles added.
    std::vector<int> Tlast;
    for (int d = 0; d < 64; d++)
    {
        Tlast.push_back(-1);
    }
    int numL = 64; // number  slots.
    int L_idx = -1;

    // keep track of covered (dead) triangles.
    std::vector<int> Xout;
    for (int d = 0; d < 64; d++)
    {
        Xout.push_back(-1);
    }
    int numX = 64;  // number of dead triangle slots.
    int X_idx = -1; // index to next dead triangle slot to be retuned to use.

    float mr = 0, mc = 0, mz = 0;
    float Mr = 0, Mc = 0, Mz = 0;

    Tri T1(0, 1, 2);
    float r0 = pts[0].r, c0 = pts[0].c, z0 = pts[0].z;
    float r1 = pts[1].r, c1 = pts[1].c, z1 = pts[1].z;
    float r2 = pts[2].r, c2 = pts[2].c, z2 = pts[2].z;

    Mr = r0 + r1 + r2;
    Mc = c0 + c1 + c2;
    Mz = z0 + z1 + z2;

    // check for colinearity
    float r01 = r1 - r0, r02 = r2 - r0;
    float c01 = c1 - c0, c02 = c2 - c0;
    float z01 = z1 - z0, z02 = z2 - z0;

    float e0 = c01 * z02 - c02 * z01;
    float e1 = -r01 * z02 + r02 * z01;
    float e2 = r01 * c02 - r02 * c01;

    if (e0 == 0 && e1 == 0 && e2 == 0)
    { // do not add a facet.
        cerr << "stop fucking me arround and give me a valid opening facet, you tit. " << endl;
        return (-1);
    }

    T1.id = 0;
    T1.er = e0;
    T1.ec = e1;
    T1.ez = e2;

    T1.ab = 1; // adjacent facet id number
    T1.ac = 1;
    T1.bc = 1;

    hull.push_back(T1);

    T1.id = 1;
    T1.er = -e0;
    T1.ec = -e1;
    T1.ez = -e2;

    T1.ab = 0;
    T1.ac = 0;
    T1.bc = 0;

    hull.push_back(T1);
    std::vector<int> xlist;
    Tri Tnew;

    for (int p = 3; p < nump; p++)
    { // add points until a non coplanar set of points is achieved.
        R3 &pt = pts[p];

        Mr += pt.r;
        mr = Mr / (p + 1);
        Mc += pt.c;
        mc = Mc / (p + 1);
        Mz += pt.z;
        mz = Mz / (p + 1);

        // find the first visible plane.
        int numh = (int)hull.size();
        int hvis = -1;
        float r = pt.r;
        float c = pt.c;
        float z = pt.z;
        xlist.clear();

        if (L_idx >= 0)
        {
            for (int L = L_idx; L >= 0; L--)
            {
                int h = Tlast[L];

                Tri &t = hull[h];
                float R1 = pts[t.a].r;
                float C1 = pts[t.a].c;
                float Z1 = pts[t.a].z;

                float dr = r - R1;
                float dc = c - C1;
                float dz = z - Z1;

                float d = dr * t.er + dc * t.ec + dz * t.ez;

                if (d > 0)
                {
                    hvis = h;
                    hull[h].keep = 0;
                    xlist.push_back(hvis);

                    //	      cerr << "y" ;

                    break;
                }
            }
        }

        if (hvis <= 0)
        {
            for (int h = numh - 1; h >= 0; h--)
            {
                Tri &t = hull[h];
                float R1 = pts[t.a].r;
                float C1 = pts[t.a].c;
                float Z1 = pts[t.a].z;

                float dr = r - R1;
                float dc = c - C1;
                float dz = z - Z1;

                float d = dr * t.er + dc * t.ec + dz * t.ez;

                if (d > 0 && hull[h].keep > 0)
                {
                    hvis = h;
                    hull[h].keep = 0;
                    xlist.push_back(hvis);

                    //cerr << "n" ;

                    break;
                }
            }
        }

        if (hvis < 0)
        {
            add_coplanar(pts, hull, p);
        }
        if (hvis >= 0)
        {
            L_idx = -1;

            // new triangular facets are formed from neighbouring invisible planes.
            int numh = (int)hull.size();
            int numx = (int)xlist.size();
            for (int x = 0; x < numx; x++)
            {
                int xid = xlist[x];
                int ab = hull[xid].ab; // facet adjacent to line ab
                Tri &tAB = hull[ab];

                float R1 = pts[tAB.a].r; // point on next triangle
                float C1 = pts[tAB.a].c;
                float Z1 = pts[tAB.a].z;

                float dr = r - R1;
                float dc = c - C1;
                float dz = z - Z1;

                float d = dr * tAB.er + dc * tAB.ec + dz * tAB.ez;

                if (d > 0)
                { // add to xlist.
                    if (hull[ab].keep == 1)
                    {
                        hull[ab].keep = 0;
                        xlist.push_back(ab);
                        numx++;
                    }
                }
                else
                { // spawn a new triangle.
                    Tnew.keep = 2;
                    Tnew.a = p;
                    Tnew.b = hull[xid].a;
                    Tnew.c = hull[xid].b;

                    Tnew.ab = -1;
                    Tnew.ac = -1;
                    Tnew.bc = ab;

                    // make normal vector.
                    float dr1 = pts[Tnew.a].r - pts[Tnew.b].r, dr2 = pts[Tnew.a].r - pts[Tnew.c].r;
                    float dc1 = pts[Tnew.a].c - pts[Tnew.b].c, dc2 = pts[Tnew.a].c - pts[Tnew.c].c;
                    float dz1 = pts[Tnew.a].z - pts[Tnew.b].z, dz2 = pts[Tnew.a].z - pts[Tnew.c].z;

                    float er = (dc1 * dz2 - dc2 * dz1);
                    float ec = -(dr1 * dz2 - dr2 * dz1);
                    float ez = (dr1 * dc2 - dr2 * dc1);

                    dr = mr - r; // points from new facet towards [mr,mc,mz]
                    dc = mc - c;
                    dz = mz - z;
                    // make it point outwards.

                    float dromadery = dr * er + dc * ec + dz * ez;

                    if (dromadery > 0)
                    {
                        Tnew.er = -er;
                        Tnew.ec = -ec;
                        Tnew.ez = -ez;
                    }
                    else
                    {
                        Tnew.er = er;
                        Tnew.ec = ec;
                        Tnew.ez = ez;
                    }

                    // try to reuse a Dead triangle.
                    int new_flag = 1, H_idx = (int)hull.size();

                    if (D_idx >= 0)
                    {
                        H_idx = Dlist[D_idx];
                        D_idx--;
                        new_flag = -1;
                    }

                    Tnew.id = H_idx;

                    // update the touching triangle tAB
                    int A = hull[xid].a, B = hull[xid].b;
                    if ((tAB.a == A && tAB.b == B) || (tAB.a == B && tAB.b == A))
                    {
                        tAB.ab = H_idx;
                    }
                    else if ((tAB.a == A && tAB.c == B) || (tAB.a == B && tAB.c == A))
                    {
                        tAB.ac = H_idx;
                    }
                    else if ((tAB.b == A && tAB.c == B) || (tAB.b == B && tAB.c == A))
                    {
                        tAB.bc = H_idx;
                    }
                    else
                    {
                        cerr << "Oh crap, the di-lithium crystals are fucked!" << endl;
                        return (-1);
                    }

                    L_idx++;
                    if (L_idx < numL)
                    {
                        Tlast[L_idx] = H_idx;
                    }
                    else
                    {
                        Tlast.push_back(H_idx);
                    }

                    if (new_flag > 0)
                    {
                        hull.push_back(Tnew);
                    }
                    else
                    {
                        hull[H_idx] = Tnew;
                    }
                }

                // second side of the struck out triangle

                int ac = hull[xid].ac; // facet adjacent to line ac
                Tri &tAC = hull[ac];

                R1 = pts[tAC.a].r; // point on next triangle
                C1 = pts[tAC.a].c;
                Z1 = pts[tAC.a].z;

                dr = r - R1;
                dc = c - C1;
                dz = z - Z1;

                d = dr * tAC.er + dc * tAC.ec + dz * tAC.ez;

                if (d > 0)
                { // add to xlist.
                    if (hull[ac].keep == 1)
                    {
                        hull[ac].keep = 0;
                        xlist.push_back(ac);
                        numx++;
                    }
                }
                else
                { // spawn a new triangle.
                    //Tnew.id = (int) hull.size();
                    Tnew.keep = 2;
                    Tnew.a = p;
                    Tnew.b = hull[xid].a;
                    Tnew.c = hull[xid].c;

                    Tnew.ab = -1;
                    Tnew.ac = -1;
                    Tnew.bc = ac;

                    // make normal vector.
                    float dr1 = pts[Tnew.a].r - pts[Tnew.b].r, dr2 = pts[Tnew.a].r - pts[Tnew.c].r;
                    float dc1 = pts[Tnew.a].c - pts[Tnew.b].c, dc2 = pts[Tnew.a].c - pts[Tnew.c].c;
                    float dz1 = pts[Tnew.a].z - pts[Tnew.b].z, dz2 = pts[Tnew.a].z - pts[Tnew.c].z;

                    float er = (dc1 * dz2 - dc2 * dz1);
                    float ec = -(dr1 * dz2 - dr2 * dz1);
                    float ez = (dr1 * dc2 - dr2 * dc1);

                    dr = mr - r; // points from new facet towards [mr,mc,mz]
                    dc = mc - c;
                    dz = mz - z;
                    // make it point outwards.

                    float dromadery = dr * er + dc * ec + dz * ez;

                    if (dromadery > 0)
                    {
                        Tnew.er = -er;
                        Tnew.ec = -ec;
                        Tnew.ez = -ez;
                    }
                    else
                    {
                        Tnew.er = er;
                        Tnew.ec = ec;
                        Tnew.ez = ez;
                    }

                    // try to reuse a Dead triangle.
                    int new_flag = 1, H_idx = (int)hull.size();

                    if (D_idx >= 0)
                    {
                        H_idx = Dlist[D_idx];
                        D_idx--;
                        new_flag = -1;
                    }

                    Tnew.id = H_idx;

                    // update the touching triangle tAC
                    int A = hull[xid].a, C = hull[xid].c;
                    if ((tAC.a == A && tAC.b == C) || (tAC.a == C && tAC.b == A))
                    {
                        tAC.ab = H_idx;
                    }
                    else if ((tAC.a == A && tAC.c == C) || (tAC.a == C && tAC.c == A))
                    {
                        tAC.ac = H_idx;
                    }
                    else if ((tAC.b == A && tAC.c == C) || (tAC.b == C && tAC.c == A))
                    {
                        tAC.bc = H_idx;
                    }
                    else
                    {
                        cerr << "Oh crap, warp drive failure, dude!" << endl;
                        return (-1);
                    }

                    L_idx++;
                    if (L_idx < numL)
                    {
                        Tlast[L_idx] = H_idx;
                    }
                    else
                    {
                        Tlast.push_back(H_idx);
                    }

                    if (new_flag > 0)
                    {
                        hull.push_back(Tnew);
                    }
                    else
                    {
                        hull[H_idx] = Tnew;
                    }

                    //hull.push_back(Tnew);
                }

                // third side of the struck out triangle

                int bc = hull[xid].bc; // facet adjacent to line ac
                Tri &tBC = hull[bc];

                R1 = pts[tBC.a].r; // point on next triangle
                C1 = pts[tBC.a].c;
                Z1 = pts[tBC.a].z;

                dr = r - R1;
                dc = c - C1;
                dz = z - Z1;

                d = dr * tBC.er + dc * tBC.ec + dz * tBC.ez;

                if (d > 0)
                { // add to xlist.
                    if (hull[bc].keep == 1)
                    {
                        hull[bc].keep = 0;
                        xlist.push_back(bc);
                        numx++;
                    }
                }
                else
                { // spawn a new triangle.
                    // Tnew.id = (int) hull.size();
                    Tnew.keep = 2;
                    Tnew.a = p;
                    Tnew.b = hull[xid].b;
                    Tnew.c = hull[xid].c;

                    Tnew.ab = -1;
                    Tnew.ac = -1;
                    Tnew.bc = bc;

                    // make normal vector.
                    float dr1 = pts[Tnew.a].r - pts[Tnew.b].r, dr2 = pts[Tnew.a].r - pts[Tnew.c].r;
                    float dc1 = pts[Tnew.a].c - pts[Tnew.b].c, dc2 = pts[Tnew.a].c - pts[Tnew.c].c;
                    float dz1 = pts[Tnew.a].z - pts[Tnew.b].z, dz2 = pts[Tnew.a].z - pts[Tnew.c].z;

                    float er = (dc1 * dz2 - dc2 * dz1);
                    float ec = -(dr1 * dz2 - dr2 * dz1);
                    float ez = (dr1 * dc2 - dr2 * dc1);

                    dr = mr - r; // points from new facet towards [mr,mc,mz]
                    dc = mc - c;
                    dz = mz - z;
                    // make it point outwards.

                    float dromadery = dr * er + dc * ec + dz * ez;

                    if (dromadery > 0)
                    {
                        Tnew.er = -er;
                        Tnew.ec = -ec;
                        Tnew.ez = -ez;
                    }
                    else
                    {
                        Tnew.er = er;
                        Tnew.ec = ec;
                        Tnew.ez = ez;
                    }

                    // try to reuse a Dead triangle.
                    int new_flag = 1, H_idx = (int)hull.size();

                    if (D_idx >= 0)
                    {
                        H_idx = Dlist[D_idx];
                        D_idx--;
                        new_flag = -1;
                    }

                    Tnew.id = H_idx;

                    // update the touching triangle tBC
                    int B = hull[xid].b, C = hull[xid].c;
                    if ((tBC.a == B && tBC.b == C) || (tBC.a == C && tBC.b == B))
                    {
                        tBC.ab = H_idx;
                    }
                    else if ((tBC.a == B && tBC.c == C) || (tBC.a == C && tBC.c == B))
                    {
                        tBC.ac = H_idx;
                    }
                    else if ((tBC.b == B && tBC.c == C) || (tBC.b == C && tBC.c == B))
                    {
                        tBC.bc = H_idx;
                    }
                    else
                    {
                        cerr << "Oh crap, rocket engine failure" << endl;
                        return (-1);
                    }

                    L_idx++;
                    if (L_idx < numL)
                    {
                        Tlast[L_idx] = H_idx;
                    }
                    else
                    {
                        Tlast.push_back(H_idx);
                    }

                    if (new_flag > 0)
                    {
                        hull.push_back(Tnew);
                    }
                    else
                    {
                        hull[H_idx] = Tnew;
                    } // hull.push_back(Tnew);
                }
            }

            numx = xlist.size();
            for (int x = 0; x < numx; x++)
            {
                // cerr << xlist[x] << " ";

                D_idx++; // keep track of all dead triangles.
                if (D_idx < numD)
                {
                    Dlist[D_idx] = xlist[x];
                }
                else
                {
                    Dlist.push_back(xlist[x]);
                    numD++;
                }
            }
            numx = 0;

            // patch up the new triangles in hull.

            int numN = (int)hull.size();
            //std::vector<Snork> norts;
            int numS = (int)norts.size();
            int nums = 0;
            Snork snort;
            //for( int q = numN-1; q>= numh; q--){

            for (int L = L_idx; L >= 0; L--)
            {
                int q = Tlast[L];

                if (hull[q].keep > 1)
                {
                    if (nums < numS)
                    {
                        norts[nums].id = q;
                        norts[nums].a = hull[q].b;
                        norts[nums].b = 1;

                        nums++;
                    }
                    else
                    {
                        snort.id = q;
                        snort.a = hull[q].b;
                        snort.b = 1;

                        norts.push_back(snort);
                        nums++;
                        numS = (int)norts.size();
                    }

                    if (nums < numS)
                    {
                        norts[nums].id = q;
                        norts[nums].a = hull[q].c;
                        norts[nums].b = 0;

                        nums++;
                    }
                    else
                    {
                        snort.a = hull[q].c;
                        snort.b = 0;
                        norts.push_back(snort);
                        nums++;
                        numS = (int)norts.size();
                    }

                    hull[q].keep = 1;
                }
            }

            sort(norts.begin(), norts.begin() + nums);
            //            int nums = (int) norts.size();

            if (nums >= 2)
            {
                for (int s = 0; s < nums - 1; s++)
                {
                    if (norts[s].a == norts[s + 1].a)
                    {
                        // link triangle sides.
                        if (norts[s].b == 1)
                        {
                            hull[norts[s].id].ab = norts[s + 1].id;
                        }
                        else
                        {
                            hull[norts[s].id].ac = norts[s + 1].id;
                        }

                        if (norts[s + 1].b == 1)
                        {
                            hull[norts[s + 1].id].ab = norts[s].id;
                        }
                        else
                        {
                            hull[norts[s + 1].id].ac = norts[s].id;
                        }
                    }
                }
            }
        }

        /*else{
	  cerr << "still in the coplanar state you fucking baboon..." << endl;
	  // rather complicated and need to add points to the 2D-hull as two faced triangles.
	  exit(0);
	  }*/

        //cerr << D_idx << " "  ;
    }

    cerr << "max triangles used " << hull.size() << endl;

    return (0);
}
