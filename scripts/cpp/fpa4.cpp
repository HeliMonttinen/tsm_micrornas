/***************************************************************************
 *   Copyright (C) 2021 by Ari Löytynoja                              *
 *   ari.loytynoja@gmail.com                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;
typedef std::vector<std::string> MyList;


class FPA3
{

    /******************** widely used variables ******************************/

    vector<int> index1;
    vector<int> index2;
    vector<int> rindex1;
    vector<int> rindex2;

    vector<int> fseq1;
    vector<int> fseq2;

    vector<int> seq1;
    vector<int> rev1;
    vector<int> seq2;

    vector<int> mask1;

    int slg;
    int fsl1;
    int fsl2;


    int sl1;
    int sl2;
    int start1;
    int end1;
    int start2;
    int end2;
    int true_start2;
    int true_end2;

    int clus_start1;
    int clus_end1;
    int clus_start2;
    int clus_end2;

    string chrom="0";
    int chrom_start = 0;

    string qry_name;
    string ref_name;

    struct Fasta_entry
    {
        string name;
        string sequence;
        int length;
    };

    struct switchPoint
    {
        int i;
        int j;
    };

    struct seqCoordinate
    {
        int pos_x;
        int pos_y;
        int matrix;
    };

    template<typename T>
    struct Array2D
    {
        private:
            int width;
            int org_width;
            int org_height;
            T * data;
        public:
            T& operator() (int x, int y) { return data[y*width + x]; }
            Array2D(const int w, const int h) : width(w), org_width(w), org_height(h) { data = new T[w*h]; }
            void resize(int nw, int nh) { if(nw*nh <= org_width*org_height) { width=nw; } else { delete [] data; data = new T[nw*nh]; org_width = nw; org_height = nh; } }
            ~Array2D() { delete [] data; }
    };

    enum Move_ptr {match=-1, xgap=-2, ygap=-3, none=-4};

    /******************** widely used variables ******************************/



    /******************** command-line argument ******************************/


    bool maximise_23_score = false;
    bool allow_23_gaps = false;
    bool maximise_score = false;
    bool maximise_length = false;
    bool verbose = false;
    bool reverse = false;
    bool swap_pair = false;
    bool scan = false;
    bool long_output = false;
    bool print_file = false;
    bool debug = false;
    bool force_overlap = false;
    bool perfect_copy = false;
    bool perfect_iupac = false;
    bool clean_rna = false;
    bool iupac = false;

    int ref_flank = 200;
    int scan_flank = 200;
    int scan_window_width = 10;
    int scan_window_limit = 2;
    int switch_flank = 200;
    int max_event_length = 1000;
    int min_length = 5;
		int gapscore = -15;

    /******************** command-line argument ******************************/





    /********************    alignment stuff    ******************************/

    void build_indeces(string *s1,string *s2)
    {

        if(s1->length() != s2->length())
        {
            cout<<"expecting aligned sequences. exiting.\n\n";
            exit(0);
        }

        if(clean_rna)
        {

            for(int i=0;i<(int)s1->length();)
            {
                if(s1->at(i) == '-' && s2->at(i) == '-')
                {
                    s1->erase(i,1);
                    s2->erase(i,1);
                }
                else
                {
                    i++;
                }
            }

            replace(s1->begin(),s1->end(),'U','T');
            replace(s1->begin(),s1->end(),'u','T');
            replace(s2->begin(),s2->end(),'U','T');
            replace(s2->begin(),s2->end(),'u','T');
        }

        slg = s1->length();

        index1.reserve(slg);
        index2.reserve(slg);
        rindex1.reserve(slg);
        rindex2.reserve(slg);

        fseq1.reserve(slg);
        fseq2.reserve(slg);

        int p1=0; int p2=0;

        for(int i=0;i<slg;i++)
        {
            index1.push_back(p1);
            index2.push_back(p2);

            if(s1->at(i) != '-')
            {
                rindex1.push_back(i);
                p1++;
            }
            if(s2->at(i) != '-')
            {
                rindex2.push_back(i);
                p2++;
            }
        }

        fsl1 = p1;
        fsl2 = p2;

        rindex1.resize(p1);
        rindex2.resize(p2);


        string alpha = "ACGTRYMKWSBDHVN";

        int ci;

        string::iterator sit = s1->begin();
        for(;sit != s1->end();sit++)
        {
            ci = alpha.find(toupper(*sit));
            if (ci>=0 && ci<15)
                fseq1.push_back(ci);
            else
                fseq1.push_back(-1);
        }

        sit = s2->begin();
        for(;sit != s2->end();sit++)
        {
            ci = alpha.find(toupper(*sit));
            if (ci>=0 && ci<15)
                fseq2.push_back(ci);
            else
                fseq2.push_back(-1);
        }

    }

    void build_sequences(string *s1,string *s2)
    {
        seq1.reserve(s1->length());
        rev1.reserve(s1->length());
        seq2.reserve(s2->length());
        mask1.reserve(s1->length());

        string alpha = "ACGTRYMKWSBDHVN";
        int revind [15] = {3,2,1,0,5,4,7,6,8,9,13,12,11,10,14};

        int ci;
        int p1=0; int p2=0;

        string::iterator sit = s1->begin();
        string::iterator sit2 = s2->begin();
        for(;sit != s1->end() && sit2 != s2->end();sit++,sit2++)
        {
            ci = alpha.find(toupper(*sit));
            if (ci>=0 && ci<15)
            {
                seq1.push_back(ci);

                int m = 0;
                if(islower(*sit))
                    m += 1;
                if(islower(*sit2))
                    m += 2;

                mask1.push_back(m);

                p1++;
            }
        }


        sit = s1->begin();
        if(reverse)
        {
            for(;sit != s1->end();sit++)
            {
                ci = alpha.find(toupper(*sit));
                if (ci>=0 && ci<15)
                    rev1.push_back(ci);

            }
        }
        else
        {
            for(;sit != s1->end();sit++)
            {
                ci = alpha.find(toupper(*sit));
                if (ci>=0 && ci<15)
                    rev1.push_back(revind[ci]);
            }
        }

        sit = s2->begin();
        for(;sit != s2->end();sit++)
        {
            ci = alpha.find(toupper(*sit));
            if (ci>=0 && ci<15)
            {
                seq2.push_back(ci);
                p2++;
            }
        }

        seq1.resize(p1);
        rev1.resize(p1);
        seq2.resize(p2);
        mask1.resize(p1);


    }


    void set_two_fragments()
    {
        start1 = 0;
        end1 = seq1.size();
        sl1 = end1-start1;


        start2 = 0;
        end2 = seq2.size();

        int flank = ref_flank;

        true_start2 = start2;
        true_end2 = end2;

        if(flank>0 && start2-flank>=0)
            start2 = start2-flank;
        if(flank>0 && end2+flank<(int)seq2.size())
            end2 = end2+flank;

        sl2 = end2-start2;

    }

    void set_alignment_fragments()
    {

        int ps = 0;
        int pe = slg-1;

        start1 = index1.at(ps);
        end1 = index1.at(pe);

        sl1 = end1-start1;

        int flank = ref_flank;

        start2 = index2.at(ps)-flank;
        end2 = index2.at(pe)+flank;

        if(start2 < 0)
            start2 = 0;

        if(end2>(int)seq2.size())
            end2 = seq2.size();

        sl2 = end2-start2;

    }

    int substitution_score(int i, int j)
    {
        if(iupac)
        {
            int smat [225] =
            {
                 10,-11,-11,-11, 10,-11, 10,-11, 10,-11,-11, 10, 10, 10,-11,
                -11, 10,-11,-11,-11, 10, 10,-11,-11, 10, 10,-11, 10, 10,-11,
                -11,-11, 10,-11, 10,-11,-11, 10,-11, 10, 10, 10,-11, 10,-11,
                -11,-11,-11, 10,-11, 10,-11, 10, 10,-11, 10, 10, 10,-11,-11,
                 10,-11, 10,-11, 10,-11, 10, 10, 10, 10, 10, 10, 10, 10,-11,
                -11, 10,-11, 10,-11, 10, 10, 10, 10, 10, 10, 10, 10, 10,-11,
                 10, 10,-11,-11, 10, 10, 10,-11, 10, 10, 10, 10, 10, 10,-11,
                -11,-11, 10, 10, 10, 10,-11, 10, 10, 10, 10, 10, 10, 10,-11,
                 10,-11,-11, 10, 10, 10, 10, 10, 10,-11, 10, 10, 10, 10,-11,
                -11, 10, 10,-11, 10, 10, 10, 10,-11, 10, 10, 10, 10, 10,-11,
                -11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,-11,
                 10,-11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,-11,
                 10, 10,-11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,-11,
                 10, 10, 10,-11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,-11,
                -11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11
            };
            if(seq1.at(i-1)>14 || seq2.at(j-1)>14)
                return -11;
            else
                return smat[seq1.at(i-1)*15+seq2.at(j-1)];
        }
        else
        {
            if(seq1.at(i-1)<4 && seq1.at(i-1) == seq2.at(j-1))
                return 10;
            else
                return -10;
        }
    }

    int char_ident(int i, int j)
    {
        if(iupac)
        {

            int smat [225] =
            {
            1,0,0,0,1,0,1,0,1,0,0,1,1,1,0,
            0,1,0,0,0,1,1,0,0,1,1,0,1,1,0,
            0,0,1,0,1,0,0,1,0,1,1,1,0,1,0,
            0,0,0,1,0,1,0,1,1,0,1,1,1,0,0,
            1,0,1,0,1,0,1,1,1,1,1,1,1,1,0,
            0,1,0,1,0,1,1,1,1,1,1,1,1,1,0,
            1,1,0,0,1,1,1,0,1,1,1,1,1,1,0,
            0,0,1,1,1,1,0,1,1,1,1,1,1,1,0,
            1,0,0,1,1,1,1,1,1,0,1,1,1,1,0,
            0,1,1,0,1,1,1,1,0,1,1,1,1,1,0,
            0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,
            1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,
            1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,
            1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
            };

            if(i<0 || j<0)
                return 0;

            if(i>15 || j>15)
                return 0;

            return smat[i*15+j];
        }
        else
        {
            return i == j;
        }
    }

    int rev_char(int i)
    {
        int revind [15] = {3,2,1,0,5,4,7,6,8,9,13,12,11,10,14};
        return revind[i];
    }

    int rev_substitution_score(int i, int j)
    {

        if(perfect_iupac)
        {
            if(iupac)
            {
                int smat [225] =
                {
                    10,-1100,-1100,-1100, 10,-1100, 10,-1100, 10,-1100,-1100, 10, 10, 10,-1100,
                   -1100, 10,-1100,-1100,-1100, 10, 10,-1100,-1100, 10, 10,-1100, 10, 10,-1100,
                   -1100,-1100, 10,-1100, 10,-1100,-1100, 10,-1100, 10, 10, 10,-1100, 10,-1100,
                   -1100,-1100,-1100, 10,-1100, 10,-1100, 10, 10,-1100, 10, 10, 10,-1100,-1100,
                    10,-1100, 10,-1100, 10,-1100, 10, 10, 10, 10, 10, 10, 10, 10,-1100,
                   -1100, 10,-1100, 10,-1100, 10, 10, 10, 10, 10, 10, 10, 10, 10,-1100,
                    10, 10,-1100,-1100, 10, 10, 10,-1100, 10, 10, 10, 10, 10, 10,-1100,
                   -1100,-1100, 10, 10, 10, 10,-1100, 10, 10, 10, 10, 10, 10, 10,-1100,
                    10,-1100,-1100, 10, 10, 10, 10, 10, 10,-1100, 10, 10, 10, 10,-1100,
                   -1100, 10, 10,-1100, 10, 10, 10, 10,-1100, 10, 10, 10, 10, 10,-1100,
                   -1100, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,-1100,
                    10,-1100, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,-1100,
                    10, 10,-1100, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,-1100,
                    10, 10, 10,-1100, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,-1100,
                   -1100,-1100,-1100,-1100,-1100,-1100,-1100,-1100,-1100,-1100,-1100,-1100,-1100,-1100,-1100
                };
                if(rev1.at(i-1)>14 || seq2.at(j-1)>14)
                    return -1100;
                else
                    return smat[rev1.at(i-1)*15+seq2.at(j-1)];

            }
            else
            {
                if(perfect_copy)
                {
                    if(rev1.at(i-1)<4 && rev1.at(i-1) == seq2.at(j-1))
                        return 10;
                    else
                        return -10000;
                }
                else
                {
                    if(rev1.at(i-1)<4 && rev1.at(i-1) == seq2.at(j-1))
                        return 10;
                    else
                        return -11;
                }
            }
        }
        else
        {
            if(iupac)
            {
                int smat [225] =
                {
                    10,-11,-11,-11, 10,-11, 10,-11, 10,-11,-11, 10, 10, 10,-11,
                   -11, 10,-11,-11,-11, 10, 10,-11,-11, 10, 10,-11, 10, 10,-11,
                   -11,-11, 10,-11, 10,-11,-11, 10,-11, 10, 10, 10,-11, 10,-11,
                   -11,-11,-11, 10,-11, 10,-11, 10, 10,-11, 10, 10, 10,-11,-11,
                    10,-11, 10,-11, 10,-11, 10, 10, 10, 10, 10, 10, 10, 10,-11,
                   -11, 10,-11, 10,-11, 10, 10, 10, 10, 10, 10, 10, 10, 10,-11,
                    10, 10,-11,-11, 10, 10, 10,-11, 10, 10, 10, 10, 10, 10,-11,
                   -11,-11, 10, 10, 10, 10,-11, 10, 10, 10, 10, 10, 10, 10,-11,
                    10,-11,-11, 10, 10, 10, 10, 10, 10,-11, 10, 10, 10, 10,-11,
                   -11, 10, 10,-11, 10, 10, 10, 10,-11, 10, 10, 10, 10, 10,-11,
                   -11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,-11,
                    10,-11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,-11,
                    10, 10,-11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,-11,
                    10, 10, 10,-11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,-11,
                   -11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11
                };
                if(rev1.at(i-1)>14 || seq2.at(j-1)>14)
                    return -11;
                else
                    return smat[rev1.at(i-1)*15+seq2.at(j-1)];

            }
            else
            {
                if(perfect_copy)
                {
                    if(rev1.at(i-1)<4 && rev1.at(i-1) == seq2.at(j-1))
                        return 10;
                    else
                        return -10000;
                }
                else
                {
                    if(rev1.at(i-1)<4 && rev1.at(i-1) == seq2.at(j-1))
                        return 10;
                    else
                        return -11;
                }
            }
        }
    }
    /********************    alignment stuff    ******************************/


    /********************    alignment output   ******************************/

    void print_switch_process(vector<seqCoordinate> *path,vector<switchPoint> *points)
    {
        string alpha = "ACGTRYMKWSBDHVN";

        string out1("");
        string out2("");
        string out2_gaps;
        string out3(" ");


        int point2 = points->at(1).j;
        int point3 = points->at(2).j;

        int p=path->size()-1;
        int site1 = path->at(p).pos_x;
        int site2 = path->at(p).pos_y;
        int site_mat = path->at(p).matrix;

        int first_site2 = site2;
        if(first_site2<0)
        {
            for(int i=p-1;i>0 && first_site2<0;i--)
                first_site2 = path->at(i).pos_y;
        }


        string qry(" ");
        string ref(" ");
        string rref(" ");

        int ref_pos = site2+1;
        int ref_end = 0;
        if(point3<site2)
        {
            qry += " ";
            ref += " ";
            rref += " ";
            out1 += " ";
            out3 += " ";

            for(int i=point3;i<site2;i++)
            {
                out1 += " ";
                out3 += " ";
                ref += alpha.at(seq2.at(i-1));
                if(seq2.at(i-1)<4)
                    rref += alpha.at(3-seq2.at(i-1));
                else
                    rref += "N";
                ref_end = i-1;
            }
        }

        out1 += "\bL ";

        while(site_mat==1)
        {
            if(site1>=0 && site1<=(int)seq1.size())
            {
                out1+=alpha.at(seq1.at(site1-1));
                qry+=alpha.at(seq1.at(site1-1));
            }
            else
            {
                out1+="-";
                qry+="-";
            }
            if(site2>=0 && site2<=(int)seq2.size())
            {
                ref+=alpha.at(seq2.at(site2-1));
                if(seq2.at(site2-1)<15)
                    rref+=alpha.at(rev_char(seq2.at(site2-1)));
                else
                    rref += "N";
                ref_pos = site2+1;
            }
            else
            {
                ref+="-";
                rref+="-";
                out2_gaps+=" ";
                out3+=" ";
            }
            ref_end = site2-1;

            p--;
            site1 = path->at(p).pos_x;
            site2 = path->at(p).pos_y;
            site_mat = path->at(p).matrix;
        }
        out1+=string(" 1");

        int first_site3 = site2;
        if(first_site3<0)
        {
            for(int i=p-1;i>0 && first_site3<0;i--)
                first_site3 = path->at(i).pos_y;
        }

        qry+="1 3";

        int last_y=0;

        while(site_mat==2)
        {
            if(site1>=0 && site1<=(int)seq1.size())
            {
                out2=alpha.at(seq1.at(site1-1))+out2;
                qry+=alpha.at(seq1.at(site1-1));
            }
						else if(site1==-1)
            {
							out2="-"+out2;
							qry+="-";
						}
            else
            {
                cout<<"error!\n";
            }

            last_y = site2;

            p--;
            site1 = path->at(p).pos_x;
            site2 = path->at(p).pos_y;
            site_mat = path->at(p).matrix;
        }


        if(last_y>1)
            out2_gaps+=" ";

        if(ref_pos<=0)
            ref_pos=1;
        for(int i=ref_pos;i<site2 && i<(int)seq2.size();i++)
        {
            ref+=alpha.at(seq2.at(i-1));
            if(seq2.at(i-1)<15)
                rref += alpha.at(rev_char(seq2.at(i-1)));
            else
                rref += "N";
            ref_end = i-1;
        }

        out2=out2_gaps+string("\b3 ")+out2+string(" 2");
        for(int i=first_site2;i<last_y-1;i++)
        {
            out2=string(" ")+out2;
        }


        for(int i=first_site2;i<site2-1;i++)
        {
            out3+=string(" ");
            if(seq2.at(i)<0)
                out3+=string(" ");
        }

        out3+=string("\b4 ");

        qry+="2 4";

        while(site_mat==3 )
        {
            if(site1>=0 && site1<=(int)seq1.size())
            {
                out3+=alpha.at(seq1.at(site1-1));
                qry+=alpha.at(seq1.at(site1-1));
            }
            else
            {
                out3+="-";
                qry+="-";
            }
            if(site2>=0 && site2<=(int)seq2.size() && site2-1>ref_end)
            {
                ref+=alpha.at(seq2.at(site2-1));
                if(seq2.at(site2-1)<15)
                    rref+=alpha.at(rev_char(seq2.at(site2-1)));
                else
                    rref += "N";
            }
            else if( site2-1>ref_end )
            {
                ref+="-";
                rref+="-";
            }

            p--;
            if(p>=0)
            {
                site1 = path->at(p).pos_x;
                site2 = path->at(p).pos_y;
                site_mat = path->at(p).matrix;
            }
            else
                break;
        }
        site2++;
        if(point2>site2)
        {
            for(;site2<=point2;site2++)
            {
                ref+=alpha.at(seq2.at(site2-1));
                if(seq2.at(site2-1)<14)
                    rref+=alpha.at(rev_char(seq2.at(site2-1)));
                else
                    rref += "N";
            }
        }
        out3+=string(" R");

        string qry0;
        for(int i=start1;i<end1;i++)
            qry0+=alpha.at(seq1.at(i));

        string ref0;
        for(int i=start2;i<end2;i++)
            ref0+=alpha.at(seq2.at(i));

        string ref00;
        for(int i=true_start2;i<true_end2;i++)
            ref00+=alpha.at(seq2.at(i));

        if(long_output && !scan)
        {
            cout<<endl<<"chr"<<chrom<<":"<<chrom_start+points->at(0).i+1<<"-"<<chrom_start+points->at(0).i-points->at(2).j+points->at(1).j+1<<endl<<endl;
        }

        if(verbose)
        {
            cout<<"Switch process:"<<"\nF1:  "<<out1<<"\nF3:  "<<out3<<"\nRF:  "<<ref<<"\nRR:  "<<rref<<"\nF2:  "<<out2<<"\n\n";
        }
    }



    void print_inversion_fragment(vector<seqCoordinate> *path,vector<switchPoint> *points,vector<seqCoordinate> *fwd_path, MyList *hits)
    {

        string alpha =    "-ACGTRYMKWSBDHVN";
        string lowalpha = "-acgtrymkwsbdhvn";
        int flanking = switch_flank;

        int start_i = max(points->at(0).i-flanking,0);
        int stop_i = points->at(0).i;

        int p=path->size()-1;

        int site1 = path->at(p).pos_x;
        int site2 = path->at(p).pos_y;

        while(site1<=start_i)
        {
            p--;
            site1 = path->at(p).pos_x;
            site2 = path->at(p).pos_y;
        }

        vector<int> seq1_frag1;
        vector<int> seq2_frag1;
        vector<bool> low_frag1;

        while(site1<=stop_i)
        {
            if(site1>0 && site1<=(int)seq1.size() &&
               site2>0 && site2<=(int)seq2.size() &&
               not char_ident(seq1.at(site1-1),seq2.at(site2-1)) )
                low_frag1.push_back(true);
            else if(site2<0)
                low_frag1.push_back(true);
            else
                low_frag1.push_back(false);

            if(site1>0 && site1<=(int)seq1.size())
                seq1_frag1.push_back(seq1.at(site1-1));
            else
                seq1_frag1.push_back(-1);

            if(site2>0 && site2<=(int)seq2.size())
                seq2_frag1.push_back(seq2.at(site2-1));
            else
                seq2_frag1.push_back(-1);

            p--;
            site1 = path->at(p).pos_x;
            site2 = path->at(p).pos_y;
        }

        vector<int> seq1_frag2;
        vector<int> seq2_frag2;
        vector<bool> low_frag2;

        stop_i = points->at(2).i;

        int sA=0; int sC=0; int sG=0; int sT=0;
        while(site1<=stop_i)
        {
            if(site1>0 && site1<=(int)seq1.size() &&
               site2>0 && site2<=(int)seq2.size() &&
              ( (reverse && not char_ident(seq1.at(site1-1), seq2.at(site2-1)) ) ||
                ( not reverse && not char_ident(seq1.at(site1-1),rev_char(seq2.at(site2-1))) ) ) )
                low_frag2.push_back(true);
            else if(site2<0)
                low_frag2.push_back(true);
            else
                low_frag2.push_back(false);

            if(site1>0 && site1<=(int)seq1.size())
            {
                seq1_frag2.push_back(seq1.at(site1-1));
                int c = seq1.at(site1-1);
                if(c==0) sA=1;
                if(c==1) sC=1;
                if(c==2) sG=1;
                if(c==3) sT=1;
            }
            else
                seq1_frag2.push_back(-1);

            if(site2>0 && site2<=(int)seq2.size())
                seq2_frag2.push_back(rev_char(seq2.at(site2-1)));
            else
                seq2_frag2.push_back(-1);

            p--;
            site1 = path->at(p).pos_x;
            site2 = path->at(p).pos_y;
        }
        int sumNuc = sA+sC+sG+sT;

        vector<int> seq1_frag3;
        vector<int> seq2_frag3;
        vector<bool> low_frag3;

        stop_i = min(points->at(3).i+flanking,(int)seq1.size());

        while(site1<stop_i && p>=0)
        {
            if(site1>0 && site1<=(int)seq1.size() &&
               site2>0 && site2<=(int)seq2.size() &&
               not char_ident(seq1.at(site1-1), seq2.at(site2-1)) )
                low_frag3.push_back(true);
            else if(site2<0)
                low_frag3.push_back(true);
            else
                low_frag3.push_back(false);

            if(site1>0 && site1<=(int)seq1.size())
                seq1_frag3.push_back(seq1.at(site1-1));
            else
                seq1_frag3.push_back(-1);

            if(site2>0 && site2<=(int)seq2.size())
                seq2_frag3.push_back(seq2.at(site2-1));
            else
                seq2_frag3.push_back(-1);

            p--;
            if(p>=0)
            {
                site1 = path->at(p).pos_x;
                site2 = path->at(p).pos_y;
            }
        }

        vector<int> seq1_fwd;
        vector<int> seq2_fwd;
        vector<bool> low_fwd;

        p=fwd_path->size()-1;

        site1 = fwd_path->at(p).pos_x;
        site2 = fwd_path->at(p).pos_y;

        while(site1<=start_i)
        {
            p--;
            site1 = fwd_path->at(p).pos_x;
            site2 = fwd_path->at(p).pos_y;
        }
        stop_i = min(points->at(3).i+flanking,(int)seq1.size());
        int epo_start1 = site1-1;
        int epo_stop1 = stop_i-1;

        while(site1<stop_i && p>=0)
        {
            if(site1>0 && site1<=(int)seq1.size() &&
               site2>0 && site2<=(int)seq2.size() &&
               not char_ident(seq1.at(site1-1), seq2.at(site2-1)) )
                low_fwd.push_back(true);
            else if(site2<0)
                low_fwd.push_back(true);
            else
                low_fwd.push_back(false);

            if(site1>0 && site1<=(int)seq1.size())
                seq1_fwd.push_back(seq1.at(site1-1));
            else
                seq1_fwd.push_back(-1);

            if(site2>0 && site2<=(int)seq2.size())
                seq2_fwd.push_back(seq2.at(site2-1));
            else
                seq2_fwd.push_back(-1);

            p--;
            if(p>=0)
            {
                site1 = fwd_path->at(p).pos_x;
                site2 = fwd_path->at(p).pos_y;
            }
        }


        float up_ident;
        float repeat_ident;
        float down_ident;
        float inv_ident;
        float fwd_ident;

        int inv_sum_length = 0;
        int inv_sum_ins = 0;
        int inv_sum_del = 0;
        int inv_sum_mis = 0;

        int sum1=0;
        for(int i=0;i<(int)seq1_frag1.size();i++)
        {
            if( char_ident(seq1_frag1.at(i),seq2_frag1.at(i)) )
                sum1++;
            else if(seq1_frag1.at(i)<0)
                inv_sum_del++;
            else if(seq2_frag1.at(i)<0)
                inv_sum_ins++;
            else
                inv_sum_mis++;

            inv_sum_length++;
         }

        if(seq1_frag1.size()>0)
            up_ident = float(sum1)/int(seq1_frag1.size());
        else
            up_ident = 0;

        int sum2=0;
        int pState=-1;
        bool hasCG=false;
        bool hasGC=false;
        int ts_mis = 0;
        int ts_del = 0;
        int ts_ins = 0;
        for(int i=0;i<(int)seq1_frag2.size();i++)
        {
            if(not reverse && char_ident(seq1_frag2.at(i),seq2_frag2.at(i)))
                sum2++;
            else if(reverse && seq1_frag2.at(i)>=0 && seq1_frag2.at(i)<14 && char_ident(seq1_frag2.at(i),rev_char(seq2_frag2.at(i))))
                sum2++;
            else if(seq1_frag2.at(i)<0)
                ts_del++;
            else if(seq2_frag2.at(i)<0)
                ts_ins++;
            else
                ts_mis++;

            inv_sum_length++;

            if(pState==1 && seq1_frag2.at(i)==2)
                hasCG=true;
            else if(pState==2 && seq1_frag2.at(i)==1)
                hasGC=true;

            pState=seq1_frag2.at(i);
        }

        inv_sum_del += ts_del;
        inv_sum_ins += ts_ins;
        inv_sum_mis += ts_mis;


        int CpG=0;
        if(hasCG)
            CpG+=1;
        if(hasGC)
            CpG+=2;

        repeat_ident = float(sum2)/int(seq1_frag2.size());


        int sum3=0;
        for(int i=0;i<(int)seq1_frag3.size();i++)
        {
            if(char_ident(seq1_frag3.at(i),seq2_frag3.at(i)))
                sum3++;
            else if(seq1_frag3.at(i)<0)
                inv_sum_del++;
            else if(seq2_frag3.at(i)<0)
                inv_sum_ins++;
            else
                inv_sum_mis++;

            inv_sum_length++;
         }

        if(seq1_frag3.size()>0)
            down_ident = float(sum3)/int(seq1_frag3.size());
        else
            down_ident = 0;

        inv_ident = float(sum1+sum2+sum3)/int(seq1_frag1.size()+seq1_frag2.size()+seq1_frag3.size());

        int sum_mis = 0;
        int sum_ins = 0;
        int sum_del = 0;

        int sum4=0;
        for(int i=0;i<(int)seq1_fwd.size();i++)
        {
            if(char_ident(seq1_fwd.at(i),seq2_fwd.at(i)))
                sum4++;
            else
            {
                if(seq1_fwd.at(i)<0)
                    sum_del++;
                else if(seq2_fwd.at(i)<0)
                    sum_ins++;
                else
                    sum_mis++;
            }
        }
        fwd_ident = float(sum4)/int(seq1_fwd.size());


        float epo_ident = 0;
        int mask_state = 0;

        int clus_ins = 0;
        int clus_del = 0;
        int clus_mis = 0;

        if(scan || print_file)
        {
            int epo_sum = 0; int epo_length = 0;
            this->fwd_compare_sequences(&epo_sum,&epo_length,epo_start1,epo_stop1);
            epo_ident = float(epo_sum)/epo_length;

            int m_start = clus_start1;
            if(m_start>0)
                m_start--;
            if(m_start>0)
                m_start--;

            if(clus_end1>=fsl1)
                clus_end1 = fsl1-1;

            int m_end = clus_end1;
            if(m_end+1<fsl1)
                m_end++;
            if(m_end+1<fsl1)
                m_end++;

            for(int i=m_start;i<m_end;i++)
            {
                if(mask1.at(i))
                {
                    if(mask1.at(i)>mask_state)
                        mask_state = mask1.at(i);
                }
            }

            int clus_start = rindex1.at(clus_start1);
            int clus_end = rindex1.at(clus_end1);

            if(clus_start1>0)
                clus_start = rindex1.at(clus_start1-1);
            if(clus_end1<(int)rindex1.size()-1)
                clus_end = rindex1.at(clus_end1+1);

            for(int i=clus_start;i<clus_end;i++)
            {
                if(fseq1.at(i)<0)
                    clus_del++;
                else if(fseq2.at(i)<0)
                    clus_ins++;
                else if(fseq1.at(i) != fseq2.at(i))
                    clus_mis++;
            }
        }

        int ts_mis_f = 0;
        int ts_del_f = 0;
        int ts_ins_f = 0;
        int ts_start = points->at(0).i;
        int ts_stop = points->at(3).i-1;



        p=fwd_path->size()-1;

        site1 = fwd_path->at(p).pos_x;
        site2 = fwd_path->at(p).pos_y;

        while(site1<=ts_start)
        {
            p--;
            site1 = fwd_path->at(p).pos_x;
            site2 = fwd_path->at(p).pos_y;
        }

        // cout<<endl;
        while(site1<=ts_stop && p>=0)
        {

            int c1=-1;
            int c2=-1;

            if(site1>0 && site1<=(int)seq1.size())
                c1 = seq1.at(site1-1);

            if(site2>0 && site2<=(int)seq2.size())
                c2 = seq2.at(site2-1);


            if(c1<0)
                ts_del_f++;
            if(c2<0)
                ts_ins_f++;

            if(c1>=0 && c2>=0 && c1!=c2)
                ts_mis_f++;

            p--;
            if(p>=0)
            {
                site1 = fwd_path->at(p).pos_x;
                site2 = fwd_path->at(p).pos_y;
            }
        }

        //

        cout<<setprecision(3);

        if(scan || print_file)
        {
            if(scan || long_output)
            {
                int c_start = rindex1.size();
                if(clus_start1<(int)rindex1.size())
                    c_start = rindex1.at(clus_start1);

                stringstream ss;
                ss <<chrom<<","<<chrom_start+clus_start1<<","<<c_start<<","<<clus_start1<<","<<clus_end1<<","
                   <<points->at(0).i<<","<<points->at(0).j<<","<<points->at(1).j<<","<<points->at(2).j<<","<<points->at(3).j<<","
                   <<up_ident<<","<<repeat_ident<<","<<down_ident<<","<<inv_ident<<","<<fwd_ident<<","<<epo_ident<<","
                   <<mask_state<<","<<sum_ins-inv_sum_ins<<","<<sum_del-inv_sum_del<<","<<sum_mis-inv_sum_mis<<","<<sumNuc<<","
                   <<CpG<<","<<clus_ins<<","<<clus_del<<","<<clus_mis<<","<<ts_ins<<","<<ts_del<<","<<ts_mis<<","<<ts_ins_f<<","<<ts_del_f<<","<<ts_mis_f<<","
                   <<inv_sum_ins<<","<<inv_sum_del<<","<<inv_sum_mis<<","<<sum_ins<<","<<sum_del<<","<<sum_mis;

                hits->push_back(ss.str());
            }
        }
        else
        {
            stringstream ss;
            ss <<points->at(0).i<<","<<points->at(0).j<<","<<points->at(1).j<<","<<points->at(2).j<<","<<points->at(3).j<<","
               <<up_ident<<","<<repeat_ident<<","<<down_ident<<","<<inv_ident<<","<<fwd_ident<<","<<sum_ins<<","<<sum_del<<","<<sum_mis<<","<<sumNuc<<","<<CpG;
            hits->push_back(ss.str());
        }


        if(verbose)
        {
            if(long_output)
            {
                cout<<"\nSwitch point 1: "<<points->at(0).j<<" ("<<points->at(0).i<<")\n";
                cout<<  "       point 2: "<<points->at(1).j<<" ("<<points->at(1).i<<")\n";
                cout<<  "       point 3: "<<points->at(2).j<<" ("<<points->at(2).i<<")\n";
                cout<<  "       point 4: "<<points->at(3).j<<" ("<<points->at(3).i<<")\n";

                cout<<"\nIdentity upstream: "<<up_ident<<endl;
                cout<<"           repeat: "<<repeat_ident<<endl;
                cout<<"       downstream: "<<down_ident<<endl;
                cout<<"        inversion: "<<inv_ident<<endl;
                cout<<"          forward: "<<fwd_ident<<endl;
            }

            stringstream qry;
            stringstream ref;
            for(int i=0;i<(int)seq1_frag1.size();i++)
            {
                if(low_frag1.at(i))
                    qry<<lowalpha.at(seq1_frag1.at(i)+1);
                else
                    qry<<alpha.at(seq1_frag1.at(i)+1);
                ref<<alpha.at(seq2_frag1.at(i)+1);
            }
            ref<<"|";
            qry<<"|";
            for(int i=0;i<(int)seq1_frag2.size();i++)
            {
                if(low_frag2.at(i))
                    qry<<lowalpha.at(seq1_frag2.at(i)+1);
                else
                    qry<<alpha.at(seq1_frag2.at(i)+1);
                ref<<alpha.at(seq2_frag2.at(i)+1);
            }
            ref<<"|";
            qry<<"|";
            for(int i=0;i<(int)seq1_frag3.size();i++)
            {
                if(low_frag3.at(i))
                    qry<<lowalpha.at(seq1_frag3.at(i)+1);
                else
                    qry<<alpha.at(seq1_frag3.at(i)+1);
                ref<<alpha.at(seq2_frag3.at(i)+1);
            }
            ref<<"\n";
            qry<<"\n";

            stringstream fwd_qry;
            stringstream fwd_ref;
            for(int i=0;i<(int)seq1_fwd.size();i++)
            {
                if(low_fwd.at(i))
                    fwd_qry<<lowalpha.at(seq1_fwd.at(i)+1);
                else
                    fwd_qry<<alpha.at(seq1_fwd.at(i)+1);
                fwd_ref<<alpha.at(seq2_fwd.at(i)+1);
            }
            fwd_ref<<"\n";
            fwd_qry<<"\n";

            if(long_output)
                cout<<"\nForward alignment:\n"<<qry_name<<" "<<fwd_qry.str()<<ref_name<<" "<<fwd_ref.str()<<endl;

            cout<<"Template-switch alignment:\n"<<qry_name<<" "<<qry.str()<<ref_name<<" "<<ref.str()<<endl;

            if(scan || print_file)
            {
                cout<<"Alignment:\n";

                int start1_epo = max(0,max( fwd_path->at(fwd_path->size()-1).pos_x-1, points->at(0).i-flanking ) );
                int ss1 = fwd_path->at(0).pos_x+1;
                if(ss1<0)ss1=(int)seq1.size();
                int stop1_epo = min( points->at(3).i+flanking,min(ss1,(int)seq1.size()-1 ));

                cout<<qry_name<<" ";
                for(int i=rindex1.at(start1_epo);i<rindex1.at(stop1_epo)-1;i++)
                    if(fseq1.at(i)==fseq2.at(i))
                       cout<<alpha.at(fseq1.at(i)+1);
                    else
                       cout<<lowalpha.at(fseq1.at(i)+1);
                cout<<endl<<ref_name<<" ";

                for(int i=rindex1.at(start1_epo);i<rindex1.at(stop1_epo)-1;i++)
                    cout<<alpha.at(fseq2.at(i)+1);
                cout<<endl<<endl;
            }
        }

    }

    /********************    alignment output   ******************************/



    /********************    alignment itself   ******************************/

    void fwd_compare_sequences(int *identical,int *length,int epo_start1,int epo_stop1)
    {
        *identical = 0;
        *length = 0;

        for(int i=rindex1.at(epo_start1);i<rindex1.at(epo_stop1);i++)
        {
            if(char_ident(fseq1.at(i), fseq2.at(i)))
                (*identical)++;
            (*length)++;
        }
    }

    void fwd_compare_sequences2(int *mis, int *del, int *ins, int *length,int epo_start1,int epo_stop1)
    {
        *mis = 0;
        *del = 0;
        *ins = 0;
        *length = 0;

        for(int i=rindex1.at(epo_start1);i<rindex1.at(epo_stop1);i++)
        {
            cout<<string("-ACGTRYMKWSBDHVN").at(fseq1.at(i)+1)<<" "<<string("-ACGTRYMKWSBDHVN").at(fseq2.at(i)+1)<<endl;

            if(char_ident(fseq1.at(i), fseq2.at(i)))
                ;
            else
            {
                if(fseq1.at(i)<0)
                    (*del)++;
                else if(fseq2.at(i)<0)
                    (*ins)++;
                else
                    (*mis)++;
            }

            (*length)++;
        }
    }

    void fwd_align_sequences(vector<seqCoordinate> *path,vector<seqCoordinate> *inv_path)
    {

        int p = inv_path->size()-1;
        int si = inv_path->at(p).pos_x;
        while(si<0 && p>=0)
            si = inv_path->at(--p).pos_x;

        p = inv_path->size()-1;
        int sj = inv_path->at(p).pos_y;
        while(sj<0 && p>=0)
            sj = inv_path->at(--p).pos_y;
        p=0;
        int ei = inv_path->at(0).pos_x;
        while(ei<0 && p<(int)inv_path->size())
            ei = inv_path->at(++p).pos_x;
        p=0;
        int ej = inv_path->at(p).pos_y;
        while(ej<0 && p<(int)inv_path->size())
            ej = inv_path->at(++p).pos_y;

        int tmp = min(si,ei);
        ei = max(si,ei);
        si = tmp;
        tmp = min(sj,ej);
        ej = max(sj,ej);
        sj = tmp;


        if(scan || print_file)
        {
            si = start1;
            ei = end1;
            sj = true_start2;
            ej = true_end2;
        }

        Array2D<int> mat1(ei-si+1,ej-sj+1);
        Array2D<int> ptr1(ei-si+1,ej-sj+1);

        int large_neg = -10000;

        for(int j=0;j<=ej-sj;j++)
        {
            // mat1
            for(int i=0;i<=ei-si;i++)
            {
                int ptr = none;
                int score = large_neg;

                if(i==0 && j==0)
                    score = 0;

                else
                {
                    if(i>0 && j>0)
                    {
                        score = mat1(i-1,j-1) + substitution_score(si+i,sj+j);
                        ptr = match;
                    }

                    if(i>0)
                    {
                        int this_score = mat1(i-1,j) + gapscore;
                        if(this_score > score)
                        {
                            score = this_score;
                            ptr = xgap;
                        }
                    }

                    if(j>0)
                    {
                        int this_score = mat1(i,j-1) + gapscore;
                        if(this_score > score)
                        {
                            score = this_score;
                            ptr = ygap;
                        }
                    }
                }

                mat1(i,j) = score;
                ptr1(i,j) = ptr;

            }
        }

        int i=ei-si;
        int j=ej-sj;

        for(;i>=0 || j>=0;)
        {
            seqCoordinate c;
            c.matrix = 1;
            c.pos_x = si+i;
            c.pos_y = sj+j;

            if(ptr1(i,j)==match)
            {
                path->push_back(c);
                i--;j--;
            }

            else if(ptr1(i,j)==xgap)
            {
                c.pos_y = -1;
                path->push_back(c);
                i--;
            }

            else if(ptr1(i,j)==ygap)
            {
                c.pos_x = -1;
                path->push_back(c);
                j--;
            }

            else if(ptr1(i,j)==none)
            {
                break;
            }
        }
    }

    void align_sequences(vector<seqCoordinate> *path,vector<switchPoint> *points,bool local=true)
    {

        int maxtermgap = 0;

        Array2D<int> mat1(sl1+1,sl2+1);
        Array2D<int> mat2(sl1+1,sl2+1);
        Array2D<int> mat3(sl1+1,sl2+1);

        Array2D<int> ptr1(sl1+1,sl2+1);
        Array2D<int> ptr2(sl1+1,sl2+1);
        Array2D<int> ptr3(sl1+1,sl2+1);

        Array2D<int> sco2(sl1+1,sl2+1);
        Array2D<int> len2(sl1+1,sl2+1);
        Array2D<int> dist2(sl1+1,sl2+1);

        Array2D<int> sco3(sl1+1,sl2+1);
        Array2D<int> len3(sl1+1,sl2+1);

        int large_neg = -10000;

        string galpha = "-ACGTRYMKWSBDHVN";


        bool debug_matrix = debug;

        if(debug_matrix)
        {
            cout<<start1<<" "<<end1<<"; "<<start2<<" "<<end2<<"; "<<true_start2<<" "<<true_end2<<" | "<<clus_start1<<" "<<clus_end1<<"; "<<clus_start2<<" "<<clus_end2<<"\n";
            for(int i=0;i<=sl1;i++)
                 cout<<" "<<galpha.at(seq1.at(start1+i)+1);
            cout<<endl;
        }


        int maxs = 0;
        int maxi = 0;
        int maxj = 0;

        for(int j=0;j<=sl2;j++)
        {
            if(debug_matrix) cout<<start2+j<<" "<<galpha.at(seq2.at(start2+j)+1);
            // mat1
            for(int i=0;i<=sl1;i++)
            {
                int ptr = none;
                int score = large_neg;

                if(i==0 && j==0)
                    score = 0;
                else if( (i==0 && j<=maxtermgap) || (j==0 && i<=maxtermgap) )
                    score = 0;
                else if(not local && force_overlap && j+start2<true_start2)
                    ;
                else if(not local && force_overlap && j+start2>clus_end2)
                    ;
                else if(not local && force_overlap && i+start1>clus_end1)
                    ;
                else if(not local &&  i==0 && j+start2==true_start2 )
                    score = 0;

                else
                {
                    if(i>0 && j>0)
                    {
                        score = mat1(i-1,j-1) + substitution_score(start1+i,start2+j);
                        ptr = match;
                    }

                    if(i>0)
                    {
                        int this_score = mat1(i-1,j) + gapscore;
                        if(this_score > score)
                        {
                            score = this_score;
                            ptr = xgap;
                        }
                    }

                    if(j>0)
                    {
                        int this_score = mat1(i,j-1) + gapscore;
                        if(this_score > score)
                        {
                            score = this_score;
                            ptr = ygap;
                        }
                    }
                }

                if(score>0 || not local)
                {
                    mat1(i,j) = score;
                    ptr1(i,j) = ptr;
                }
                else if(local)
                {
                    mat1(i,j) = 0;
                    ptr1(i,j) = none;
                }

                if(score>maxs)
                {
                    maxs=score;
                    maxi=i;
                    maxj=j;
                }

                if(debug_matrix) cout<<" "<<score;
            }
            if(debug_matrix) cout<<"\n";
        }

        if(debug_matrix) cout<<"\n";
        if(debug_matrix) cout<<maxs<<" "<<maxi<<" "<<maxj<<"\n";

        for(int i=0;i<=sl1;i++  )
        {

            int pscore = large_neg;
            int pk = -1;

            // from mat1
            for(int k=0;k<=sl2;k++)
            {
             	if(mat1(i-1,k) > pscore)
               	{
                     pscore = mat1(i-1,k);
                     pk = k;
              	}
            }


	          // mat2
            for(int j=sl2;j>=0;j--)
            {

                int score = large_neg;
                int subst = 0;
                int ptr = none;
                int len = 0;
                int sco = 0;
                int dist = sl2;

                if(i>0 && j>0)
                {
                    ptr = match;
                    subst = rev_substitution_score(start1+i,start2+j);

                    if(maximise_23_score)
                    {
                      if(j<sl2 && sco2(i-1,j+1)+subst>0)
                      {
                          score = subst + mat2(i-1,j+1);
                          len = len2(i-1,j+1)+1;
                          sco = sco2(i-1,j+1)+subst;
                          dist = dist2(i-1,j+1);
                      }

                      if(allow_23_gaps)
                      {

												// if(j<sl2 && sco2(i-1,j)+gapscore>sco)
												if(j<sl2 && mat2(i-1,j)+gapscore>score)
                        {
                            score = mat2(i-1,j)+gapscore;
                            len = len2(i-1,j)+1;
                            sco = sco2(i-1,j)+gapscore;
                            dist = dist2(i-1,j);
                            ptr = xgap;
                        }

												// if(j<sl2 && sco2(i,j+1)+gapscore>sco)
												if(j<sl2 && mat2(i,j+1)+gapscore>score)
                        {
                            score = mat2(i,j+1)+gapscore;
                            len = len2(i,j+1)+1;
                            sco = sco2(i,j+1)+gapscore;
                            dist = dist2(i,j+1);
                            ptr = ygap;
                        }


                      }
                    }
                    else
                    {
                        if(j<sl2)
                        {
                            score = subst + mat2(i-1,j+1);
                            len = len2(i-1,j+1)+1;
                            sco = sco2(i-1,j+1)+subst;
                            dist = dist2(i-1,j+1);
                        }
                    }

                    // from mat1
                    if(maximise_23_score)
                    {
                        if(subst>0 && subst + pscore > score)
                        //if(subst + pscore > score)
                        {
                            score = subst + pscore;
                            ptr = pk;
                            len = 1;
                            sco = subst;
                            dist = abs(j-pk);
                        }
                    }
                    else
                    {
                        if(subst + pscore > score)
                        {
                           	score = subst + pscore;
                           	ptr = pk;
                            len = 1;
                            sco = subst;
                            dist = abs(j-pk);
                        }
                    }
                }
                mat2(i,j) = score;
                ptr2(i,j) = ptr;

                len2(i,j) = len;
                sco2(i,j) = sco;

                dist2(i,j) = dist;


								//if(debug_matrix) cout<<" "<<score<<";"<<ptr;

								if(debug_matrix) {
									if(score<0)
										cout<<" -;"<<ptr;
									else
										cout<<" "<<score<<";"<<ptr;
								}
            }
            if(debug_matrix) cout<<"\n";
        }
        if(debug_matrix) cout<<"\n";

        int max23 = large_neg;
        int max23_i = -1;
        // int max23_j = -1;

        for(int j=0;j<=sl2;j++)
        {
            for(int i=0;i<=sl1;i++)
            {
                if(sco2(i,j)>max23)
                {
                  max23 = sco2(i,j);
                  max23_i = i;
                  // max23_j = j;
                }
            }
        }

        /*
        for(int j=0;j<=sl2;j++)
        {
            for(int i=0;i<=sl1;i++)
            {
                if(sco2(i,j)>0)
                  cout<<" "<<sco2(i,j);
                else
                  cout<<" .";
            }
            cout<<endl;
        }
        cout<<endl;
        */
        // cout<<max23<<" "<<max23_i<<" "<<max23_j<<" "<<endl;

        maxs = 0;
        maxi = 0;
        maxj = 0;

        for(int i=0;i<=sl1;i++)
        {

          for(int j=0;j<=sl2;j++)
            {
                if(debug_matrix) cout<<start2+j<<" "<<galpha.at(seq2.at(start2+j)+1);
                // mat3
                int score = large_neg;
                int ptr = none;
                int len = 0;
                int sco = 0;

                if(i>0 && j>0)
                {
                    if(not local && force_overlap && j+start2<clus_start2)
                        ;
                    else if(not local && force_overlap && i+start1<clus_start1)
                        ;
                    else if(not local && force_overlap && j+start2>true_end2)
                        ;
                    else if(maximise_23_score && i<max23_i+1)
                        ;
                    else
                    {
                        ptr = match;
                        int subst = substitution_score(start1+i,start2+j);
                        score = subst + mat3(i-1,j-1);

                        len = len3(i-1,j-1);
                        sco = sco3(i-1,j-1);

                        int this_score = mat3(i-1,j) + gapscore;
                        if(this_score > score)
                        {
                            score = this_score;
                            ptr = xgap;

                            len = len3(i-1,j);
                            sco = sco3(i-1,j);
                        }

                        this_score = mat3(i,j-1) + gapscore;
                        if(this_score > score)
                        {
                            score = this_score;
                            ptr = ygap;

                            len = len3(i,j-1);
                            sco = sco3(i,j-1);
                        }
                        if(i==sl1 && j>=true_end2)
                        {
                            this_score = mat3(i,j-1);
                            if(this_score > score)
                            {
                                score = this_score;
                                ptr = ygap;

                                len = len3(i,j-1);
                                sco = sco3(i,j-1);
                            }
                        }

                        // from mat2
                        int dist = sl2;
                        if(maximise_23_score)
                        {
                          if(i==max23_i+1)
                          {
                              for(int k=0;k<=sl2;k++)
                              {
                                  //
                                  if(subst + mat2(i-1,k) > score ||
                                     ( subst + mat2(i-1,k) == score && dist2(i-1,k)<dist) )
                                    {
                                        score = subst + mat2(i-1,k);
                                        ptr = k;

                                        len = len2(i-1,k);
                                        sco = sco2(i-1,k);
                                        dist = dist2(i-1,k);

                                        if(score>maxs)
                                        {
                                            maxs=score;
                                            maxi=i;
                                            maxj=k;
                                        }

                                    }
                                }
                            }
                        }
                        else
                        {
                            for(int k=0;k<=sl2;k++)
                            {
                                if(subst + mat2(i-1,k) > score ||
                                   ( subst + mat2(i-1,k) == score && dist2(i-1,k)<dist) )
                                {
                                    score = subst + mat2(i-1,k);
                                    ptr = k;

                                    len = len2(i-1,k);
                                    sco = sco2(i-1,k);
                                    dist = dist2(i-1,k);

                                    if(score>maxs)
                                    {
                                        maxs=score;
                                        maxi=i;
                                        maxj=k;
                                    }

                                }
                            }
                        }
                    }
                }
                if(score>0 || not local)
                {
                    mat3(i,j) = score;
                    ptr3(i,j) = ptr;

                    len3(i,j) = len;
                    sco3(i,j) = sco;
                }
                else if(local)
                {
                    mat3(i,j) = 0;
                    ptr3(i,j) = none;

                    len3(i,j) = 0;
                    sco3(i,j) = 0;
                }
            if(debug_matrix) cout<<" "<<score;
            }
            if(debug_matrix) cout<<"\n";
        }
        if(debug_matrix) cout<<"\n";

        if(debug_matrix) cout<<maxs<<" "<<maxi<<" "<<maxj<<"\n";

        int max_i=-1;
        int max_j=-1;
        int max_end = large_neg;
        int max_len = large_neg;
        int max_sco = large_neg;

        for(int i=sl1;i>0;i--)
        {
            for(int j=sl2;j>0;j--)
            {
                if(maximise_score)
                {
                    if(sco3(i,j)>max_sco)
                    {
                        max_i = i;
                        max_j = j;
                        max_end = mat3(i,j);
                        max_len = len3(i,j);
                        max_sco = sco3(i,j);
                    }
                }
                else if (maximise_length)
                {
                    if(len3(i,j)>max_len)
                    {
                        max_i = i;
                        max_j = j;
                        max_end = mat3(i,j);
                        max_len = len3(i,j);
                        max_sco = sco3(i,j);
                    }
                }
                else
                {
                    if(mat3(i,j)>max_end)
                    {
                        max_i = i;
                        max_j = j;
                        max_end = mat3(i,j);
                        max_len = len3(i,j);
                        max_sco = sco3(i,j);
                    }
                }
            }
        }

        int i=max_i;
        int j=max_j;
        bool ptr_found = false;

        // cout<<i<<" "<<j<<endl;

        for(;i>=0 || j>=0;)
        {
            seqCoordinate c;
            c.matrix = 3;
            c.pos_x = start1+i;
            c.pos_y = start2+j;

            if(ptr3(i,j)==match)
            {
                path->push_back(c);

                i--;j--;
            }

            else if(ptr3(i,j)==xgap)
            {
                c.pos_y = -1;
                path->push_back(c);

                i--;
            }

            else if(ptr3(i,j)==ygap)
            {
                c.pos_x = -1;
                path->push_back(c);

                j--;
            }

            else if(ptr3(i,j)==none)
            {
                cout<<"ptr3 = none. weird. exiting.\n";
                return;
            }

            else
            {
                path->push_back(c);

                /*here -> correct*/
                points->at(3).i = start1+i;
                points->at(3).j = start2+j;

                j = ptr3(i,j);
                i--;

                points->at(2).i = start1+i;
                points->at(2).j = start2+j;

                ptr_found = true;
                break;
            }
        }

        if(!ptr_found)
        {
            cout<<"backtracking 1st path failed. exiting.\n";
	    return;
        }

        ptr_found = false;

        for(;i>=0 || j<=sl2;)
        {
            seqCoordinate c;
            c.matrix = 2;
            c.pos_x = start1+i;
            c.pos_y = start2+j;

            if(ptr2(i,j)==match)
            {
                path->push_back(c);
                i--;j++;
            }

            else if(ptr2(i,j)==xgap)
            {
                c.pos_y = -1;
                path->push_back(c);
                i--;
            }

            else if(ptr2(i,j)==ygap)
            {
                c.pos_x = -1;
                path->push_back(c);
                j++;
            }

            else if(ptr2(i,j)==none)
            {
                cout<<"ptr2 = none. weird. exiting.\n";
                return;
            }

            else
            {
								path->push_back(c);

                points->at(1).i = start1+i;
                points->at(1).j = start2+j;

                j = ptr2(i,j);
                i--;

                points->at(0).i = start1+i;
                points->at(0).j = start2+j;

                ptr_found = true;
                break;
            }
        }

        if(!ptr_found)
        {
            cout<<"backtracking 2nd path failed. exiting.\n";
            return;
        }

        for(;i>=0 || j>=0;)
        {
            seqCoordinate c;
            c.matrix = 1;
            c.pos_x = start1+i;
            c.pos_y = start2+j;

            if(ptr1(i,j)==match)
            {
                path->push_back(c);
                i--;j--;
            }

            else if(ptr1(i,j)==xgap)
            {
                c.pos_y = -1;
                path->push_back(c);
                i--;
            }

            else if(ptr1(i,j)==ygap)
            {
                c.pos_x = -1;
                path->push_back(c);
                j--;
            }

            else if(ptr1(i,j)==none)
            {
                break;
            }
        }
    }

    /********************    alignment itself   ******************************/



    /********************     alignment scan    ******************************/

    void scan_alignment(vector<seqCoordinate> *path,vector<switchPoint> *points,vector<seqCoordinate> *fwd_path, MyList *hits)
    {
        int width = scan_window_width;
        int limit = scan_window_limit;
        int flank = scan_flank;

        int sum = 0;
        int ss = 0;
        int se = slg;


        int max_length = max_event_length;
        vector<pair<int,int> > clusters;


        int i = ss;

        int c_first=-1, c_start=-1, c_end=-1;

        for(;i<width;i++)
        {
            if(not char_ident(fseq1.at(i), fseq2.at(i)))
            {
                sum++;

                if(sum==1)
                {
                    c_first = i;
                }
            }


            if(sum==limit)
            {
                c_start = i-limit+2;
                if(c_first < c_start)
                    c_start = c_first;

            }
        }

        for(;i<se;i++)
        {

            if((not char_ident(fseq1.at(i-width),fseq2.at(i-width))) && sum>0)
                sum--;
            if(not char_ident(fseq1.at(i), fseq2.at(i)))
            {
                sum++;
                if(sum==1)
                {
                    c_first = i;
                }
            }

            if(sum>=limit)
            {
                int start_pos1 = index1.at(0);
                int start_pos2 = index2.at(0);
                if(i-flank-limit>0)
                {
                    start_pos1 = index1.at(i-flank-limit);
                    start_pos2 = index2.at(i-flank-limit);
                }

                int end_pos1 = index1.at(i);
                int end_pos2 = index2.at(i);


                c_start = i-limit+1;
                if(c_first < c_start)
                    c_start = c_first;

                clus_start1 = index1.at(c_start);
                clus_start2 = index2.at(c_start);

                if(c_start>0)
                {
                    clus_start1 = index1.at(c_start-1);
                    clus_start2 = index2.at(c_start-1);
                }

                i++;
                for(;i<se;i++)
                {
                    if( (not char_ident( fseq1.at(i-width), fseq2.at(i-width))) && sum>0)
                        sum--;
                    if(not char_ident(fseq1.at(i), fseq2.at(i)))
                        sum++;

                    if(i-width+flank<slg)
                        end_pos1 = index1.at(i-width+flank);
                    else
                        end_pos1 = index1.at(slg-1);

                    if(i-width+flank<slg)
                       end_pos2 = index2.at(i-width+flank);
                    else
                        end_pos2 = index2.at(slg-1);

                    if(sum==0)
                    {

                        c_end = i-width+1;

                        clus_end1 = index1.at(c_end);
                        clus_end2 = index2.at(c_end);

                        break;
                    }
                }
                if(sum>0)
                {

                    c_end = i-width+1;

                    clus_end1 = index1.at(c_end);
                    clus_end2 = index2.at(c_end);
                }

                if(clus_end1-clus_start1<=max_length && clus_end2-clus_start2<=max_length && clus_start1<=clus_end1)
                {
                    clusters.push_back(make_pair(c_start,c_end));
                    this->check_scan_position(path,points,fwd_path,start_pos1,end_pos1,start_pos2,end_pos2,hits);
                }
            }
        }
    }



    void check_scan_position(vector<seqCoordinate> *path,vector<switchPoint> *points,vector<seqCoordinate> *fwd_path,int min1,int max1,int min2,int max2, MyList *hits)
    {
        int max_length = max_event_length;

        start1 = min1;
        start2 = min2;
        end1 = max1;
        end2 = max2;
        true_start2 = min2;
        true_end2 = max2;
        if(end1-start1>max_length || end2-start2>max_length)
            return;

        int flank = ref_flank;

        if(flank>0 && start2-flank>=0)
            start2 = start2-flank;
        if(flank>0 && end2+flank<(int)seq2.size())
            end2 = end2+flank;

        sl1 = end1-start1;
        sl2 = end2-start2;

        path->clear();
        fwd_path->clear();
        this->align_sequences(path,points,false);
        this->fwd_align_sequences(fwd_path,path);


        if(points->at(1).j-points->at(2).j+1>=min_length)
        {
            if(!long_output)
            {
                this->print_switch_process(path,points);
                this->print_inversion_fragment(path,points,fwd_path,hits);
            }
            else
            {
                this->print_inversion_fragment(path,points,fwd_path,hits);
                this->print_switch_process(path,points);
            }
        }
    }
    /********************     alignment scan    ******************************/




    /********************     print events      ******************************/

    void print_events(vector<seqCoordinate> *path,vector<switchPoint> *points,vector<seqCoordinate> *fwd_path,string filename)
    {
        verbose = true;
        int flank = scan_flank;

        ifstream input(filename.c_str(), ios::in);
        if (!input) { return; }

        string temp;
        getline(input, temp, '\n');
        bool has_line = find_if(temp.begin(), temp.end(), [](char c) { return isalpha(c); }) == temp.end();

        while(!input.eof())
        {
            if(!has_line)
                getline(input, temp, '\n');

            if(temp.size()==0)
                break;

            has_line = false;

            string chrom = temp.substr(0,temp.find(','));
            temp = temp.substr(temp.find(',')+1);

            if(chrom=="scan finished")
                break;

            stringstream tempstr(temp);

            int chrom_site,aligned,clus_st1,clus_en1,sp1_qry,sp1_ref,sp2_ref,sp3_ref,sp4_ref;
            char sep;
            tempstr>>chrom_site>>sep>>aligned>>sep>>clus_st1>>sep>>clus_en1>>sep>>sp1_qry>>sep>>sp1_ref>>sep>>sp2_ref>>sep>>sp3_ref>>sep>>sp4_ref>>sep>>temp;

            if(long_output)
                cout<<chrom<<","<<tempstr.str()<<endl;

            int start_pos1 = index1.at(0);
            int start_pos2 = index2.at(0);
            if(clus_st1>flank)
            {
                start_pos1 = index1.at(rindex1.at(clus_st1)-flank);
                start_pos2 = index2.at(rindex1.at(clus_st1)-flank);
            }

            int end_pos1 = index1.at(rindex1.at(clus_en1));
            int end_pos2 = index2.at(rindex1.at(clus_en1));
            if(rindex1.at(clus_en1)+flank<(int)index1.size())
                end_pos1 = index1.at(rindex1.at(clus_en1)+flank);
            else
                end_pos1 = index1.at(index1.size()-1);

            if(rindex1.at(clus_en1)+flank<(int)index2.size())
                end_pos2 = index2.at(rindex1.at(clus_en1)+flank);
            else
                end_pos2 = index2.at(index2.size()-1);

            clus_start1 = clus_st1;
            clus_start2 = index2.at(rindex1.at(clus_st1));

            clus_end1 = clus_en1;
            clus_end2 = index2.at(rindex1.at(clus_en1));

            MyList hits;
            this->check_scan_position(path,points,fwd_path,start_pos1,end_pos1,start_pos2,end_pos2,&hits);

        }
    }


   /********************     print events      ******************************/

    void delete_arrays(){
        index1.clear();
        index2.clear();
        rindex1.clear();
        rindex2.clear();

        fseq1.clear();
        fseq2.clear();

        seq1.clear();
        seq2.clear();
        rev1.clear();
    }

public:
    FPA3() {}
    FPA3(const string &name) {}


    MyList scan_two(string rseq, string qseq, bool verb=false)
    {

        if(index1.size()>0)
            delete_arrays();

        MyList hits;
        scan = true;

        Fasta_entry ref;
        Fasta_entry qry;

        verbose = verb;

        qry_name = "qry";
        ref_name = "ref";

        ref.sequence = rseq;
        ref.name = "ref";

        qry.sequence = qseq;
        qry.name = "qry";

        if(ref.sequence.length()==0 || qry.sequence.length()==0)
        {
            cout<<"Error in sequence input. Exiting."<<endl;
            return hits;
        }

        vector<seqCoordinate> path;
        vector<seqCoordinate> fwd_path;
        vector<switchPoint> points;
        for(int i=0;i<4;i++)
        {
            switchPoint sp;
            points.push_back(sp);
        }

        if(verbose)
        {
            cout<<"\nQuery:     "<<qry.name<<endl;
            cout<<"Reference: "<<ref.name<<endl<<endl;
        }

        this->build_indeces(&qry.sequence,&ref.sequence);
        this->build_sequences(&qry.sequence,&ref.sequence);

        this->scan_alignment(&path,&points,&fwd_path,&hits);

        return hits;
    }

    void print_hit(string rseq, string qseq, string temp, bool verb=false)
    {

        if(index1.size()>0)
            delete_arrays();

        verbose = verb;
        int flank = scan_flank;

        Fasta_entry ref;
        Fasta_entry qry;

        verbose = verb;

        qry_name = "qry";
        ref_name = "ref";

        ref.sequence = rseq;
        ref.name = "ref";

        qry.sequence = qseq;
        qry.name = "qry";

        if(swap_pair)
        {
            ref.sequence = qseq;
            qry.sequence = rseq;
        }

        if(ref.sequence.length()==0 || qry.sequence.length()==0)
        {
            cout<<"Error in sequence input. Exiting."<<endl;
            return;
        }


        vector<seqCoordinate> path;
        vector<seqCoordinate> fwd_path;
        vector<switchPoint> points;
        for(int i=0;i<4;i++)
        {
            switchPoint sp;
            points.push_back(sp);
        }

        this->build_indeces(&qry.sequence,&ref.sequence);
        this->build_sequences(&qry.sequence,&ref.sequence);

        string chrom = temp.substr(0,temp.find(','));
        temp = temp.substr(temp.find(',')+1);

        stringstream tempstr(temp);

        int chrom_site,aligned,clus_st1,clus_en1,sp1_qry,sp1_ref,sp2_ref,sp3_ref,sp4_ref;
        char sep;
        tempstr>>chrom_site>>sep>>aligned>>sep>>clus_st1>>sep>>clus_en1>>sep>>sp1_qry>>sep>>sp1_ref>>sep>>sp2_ref>>sep>>sp3_ref>>sep>>sp4_ref>>sep>>temp;

        if(long_output)
            cout<<chrom<<","<<tempstr.str()<<endl;

        int start_pos1 = index1.at(0);
        int start_pos2 = index2.at(0);
        if(clus_st1>flank)
        {
            start_pos1 = index1.at(rindex1.at(clus_st1)-flank);
            start_pos2 = index2.at(rindex1.at(clus_st1)-flank);
        } else {
            start_pos1 = index1.at(rindex1.at(1));
            start_pos2 = index2.at(rindex1.at(1));
        }

        int end_pos1 = index1.at(rindex1.at(clus_en1));
        int end_pos2 = index2.at(rindex1.at(clus_en1));
        if(rindex1.at(clus_en1)+flank<(int)index1.size())
            end_pos1 = index1.at(rindex1.at(clus_en1)+flank);
        else
            end_pos1 = index1.at(index1.size()-1);

        if(rindex1.at(clus_en1)+flank<(int)index2.size())
            end_pos2 = index2.at(rindex1.at(clus_en1)+flank);
        else
            end_pos2 = index2.at(index2.size()-1);

        clus_start1 = clus_st1;
        clus_start2 = index2.at(rindex1.at(clus_st1));

        clus_end1 = clus_en1;
        clus_end2 = index2.at(rindex1.at(clus_en1));

        print_file = true;

        MyList hits;
        this->check_scan_position(&path,&points,&fwd_path,start_pos1,end_pos1,start_pos2,end_pos2,&hits);
    }


    void set_bool(string par, bool val)
    {
        if(par=="maximise_23_score")
            maximise_23_score = val;

        if(par=="allow_23_gaps")
              allow_23_gaps = val;

        if(par=="maximise_score")
            maximise_score = val;

        if(par=="maximise_length")
            maximise_length = val;

        if(par=="verbose")
            verbose = val;

        if(par=="reverse")
            reverse = val;

        if(par=="swap_pair")
            swap_pair = val;

        if(par=="long_output")
            long_output = val;

        if(par=="force_overlap")
            force_overlap = val;

        if(par=="debug")
            debug = val;

        if(par=="perfect_copy")
            perfect_copy = val;

        if(par=="perfect_iupac")
            perfect_iupac = val;

        if(par=="clean_rna")
            clean_rna = val;

        if(par=="iupac")
            iupac = val;
    }

    void set_int(string par, int val)
    {
        if(par=="ref_flank")
            ref_flank = val;

        if(par=="scan_flank")
            scan_flank = val;

        if(par=="scan_window_width")
            scan_window_width = val;

        if(par=="scan_window_limit")
            scan_window_limit = val;

        if(par=="switch_flank")
            switch_flank = val;

        if(par=="max_event_length")
            max_event_length = val;

				if(par=="min_length")
            min_length = val;

				if(par=="gapscore")
            gapscore = val;
    }



};


#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(MyList);

PYBIND11_MODULE(fpa_ext62,m)
{

    py::bind_vector<MyList>(m, "hits");

    py::class_<FPA3>(m,"FPA2")
            .def(py::init())
            .def(py::init<const std::string &>())
            .def("scan_two", &FPA3::scan_two)
            .def("print_hit", &FPA3::print_hit)
            .def("set_int", &FPA3::set_int)
            .def("set_bool", &FPA3::set_bool)
    ;
}


/*
 * c++ -O3 -Wall -shared -std=c++11 -fPIC -I/usr/include/python3.6m -I/usr/local/include/python3.6 \
 * fpa2.cpp -o fpa_ext2.cpython-36m-x86_64-linux-gnu.so
 */

/*
 * import fpa_ext2
 * s1="ACGAGATCGATCCTCTCAGAGAGCTGACCGAGATGACGA"
 * s2="ACGAGATCGATCCTCTCTGAGAGCTGACCGAGATGACGA"
 *
 * f = fpa_ext2.FPA2()
 *
 * f.set_bool("force_overlap",True)
 * f.set_bool("iupac",True)
 * f.set_int("min_length",5)
 * f.set_int("scan_window_limit",1)
 *
 * hits1 = f.scan_two(s1,s2,False)
 *
 * for h in hits1:
 * 	f.print_hit(s1,s2,h,True)
 *
 */
