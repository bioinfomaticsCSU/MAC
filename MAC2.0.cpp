# include<bits/stdc++.h>
# include<sys/stat.h>
using namespace std;
# define bug puts("H");
# define SZ(x) (int)x.size()
# define pb push_back
# define mp make_pair
# define fi first
# define se second
# define rep(i,a,n) for (int i=a; i<=n; ++i)
# define per(i,a,n) for (int i=a; i>=n; --i)
# define all(x) (x).begin(), (x).end()
# define INF 1000000000
typedef vector<int> Vi;
typedef pair<int, int> Pi;

const string input_folder = "./input", output_folder = "./output", temp_folder = "./temp";
string qcontig_file, rcontig_file;
struct Gene{int l1, r1, l2, r2, id;};
vector<Gene> gene;
vector<Vi> qcontig, rcontig;
Vi qpre, qnxt, rpre, rnxt, qds, rds;
vector<bool> qtelo, rtelo, vis;
vector<string> contig;
int q_size = -1, r_size = -1, gene_num = 0;

void checkFile(int argc, char * argv[]) {
    if (argc != 3) {
        puts("error: Wrong number of parameters!");
        exit(0);
    }
    qcontig_file = argv[1]; rcontig_file = argv[2];
    int qlen = SZ(qcontig_file), rlen = SZ(rcontig_file);
    if (qlen <= 3 || rlen <= 3 || qcontig_file.substr(qlen-3, 3) != ".fa" || rcontig_file.substr(rlen-3, 3) != ".fa") {
        puts("error: Parameters format error!");
        exit(0);
    }
    qcontig_file = input_folder + "/" + qcontig_file; rcontig_file = input_folder + "/" + rcontig_file;
    struct stat buffer;
    if (stat(qcontig_file.c_str(), &buffer) || stat(rcontig_file.c_str(), &buffer)) {
        puts("error: Contig File does not exist!");
        exit(0);
    }
}

void cleanContig() {
    ifstream in1(qcontig_file.c_str()), in2(rcontig_file.c_str());
    qcontig_file += ".fasta"; rcontig_file += ".fasta";
    ofstream out1(qcontig_file.c_str()), out2(rcontig_file.c_str());
    int num = 0;
    string line;
    while (getline(in1, line)) {
        if (line[0] == '>') {
            if (num) out1 << endl;
            out1 << ">" << ++num << endl;
        }
        else out1 << line;
    }
    num = 0; in1.close(); out1.close();
    while (getline(in2, line)) {
        if (line[0] == '>') {
            if (num) out2 << endl;
            out2 << ">" << ++num << endl;
        }
        else out2 << line;
    }
    in2.close(); out2.close();
}

bool comp1(int a, int b) {return gene[abs(a)-1].l1 < gene[abs(b)-1].l1;}
bool comp2(int a, int b) {return gene[abs(a)-1].l2 < gene[abs(b)-1].l2;}

Vi findNumbers(string & line, Vi & inc) {
    Vi res, numarr;
    int sz = SZ(line), num = 0;
    rep(j,0,sz-1) {
        if (isdigit(line[j])) {
            num = num*10+(line[j]-'0');
            if (j+1 == sz || !isdigit(line[j+1])) numarr.pb(num), num = 0;
        }
    }
    for (auto i: inc) res.pb(numarr[i-1]);
    return res;
}

void addGene(int contig_id, int gene_id, int setid) {
    if (setid == 1) {
        rep(i,q_size+1,contig_id) qcontig.pb(Vi());
        q_size = max(q_size, contig_id);
        qcontig[contig_id].pb(gene_id);
    }
    else {
        rep(i,r_size+1,contig_id) rcontig.pb(Vi());
        r_size = max(r_size, contig_id);
        rcontig[contig_id].pb(gene_id);
    }
}

void buildContig() {
    string coords_path = temp_folder + "/cssseq.coords";
    ifstream coords(coords_path.c_str());
    string line;
    Vi inc = {1,2,3,4,13,14};
    int invaild_line = 5;
    while (invaild_line && getline(coords, line)) --invaild_line;
    while (getline(coords, line)) {
        Vi align_info = findNumbers(line, inc);
        gene.pb(Gene{min(align_info[2],align_info[3]), max(align_info[2],align_info[3]), min(align_info[0],align_info[1]), max(align_info[0],align_info[1]), align_info[5]-1});
        ++gene_num;
        addGene(align_info[4], gene_num*(align_info[2]>align_info[3]?-1:1), 2);
        addGene(align_info[5], gene_num, 1);
    }
    coords.close();
    rep(i,0,q_size) sort(all(qcontig[i]), comp1);
    rep(i,0,r_size) sort(all(rcontig[i]), comp2);
    rep(i,0,q_size) {
        int sz = SZ(qcontig[i]);
        rep(j,0,sz-1) {
            if (j == 0) gene[qcontig[i][j]-1].l1 = min(gene[qcontig[i][j]-1].l1, 1);
            else gene[qcontig[i][j]-1].l1 = min(gene[qcontig[i][j]-1].l1, gene[qcontig[i][j-1]-1].r1+1);
            if (j == sz-1) gene[qcontig[i][j]-1].r1 = INF;
        }
    }
}

int Mapping(int x){return x < 0 ? gene_num - x : x;}
int AntiMapping(int x){return x > gene_num ? gene_num-x : x;}
int Flipping(int x, int sz){return x >= sz ? 2*sz-1-x : x;}
int Sgn(int x, int sz){return x < sz ? 1 : -1;}

void initLS() {
    qpre.resize(2*gene_num+1); rpre.resize(2*gene_num+1); rnxt.resize(2*gene_num+1); qnxt.resize(2*gene_num+1); qds.resize(2*gene_num+1); rds.resize(2*gene_num+1); qtelo.resize(2*gene_num+1,false); rtelo.resize(2*gene_num+1,false);
    rep(i,0,q_size) {
        int sz = SZ(qcontig[i]);
        if (sz == 0) continue;
        rep(j,0,2*sz-1) {
            if (j == 0) qpre[Mapping(qcontig[i][0])] = Mapping(-qcontig[i][0]), qnxt[Mapping(-qcontig[i][0])] = Mapping(qcontig[i][0]);
            else qpre[Mapping(qcontig[i][Flipping(j,sz)]*Sgn(j,sz))] = Mapping(qcontig[i][Flipping(j-1,sz)]*Sgn(j-1,sz)), qnxt[Mapping(qcontig[i][Flipping(j-1,sz)]*Sgn(j-1,sz))] = Mapping(qcontig[i][Flipping(j,sz)]*Sgn(j,sz));
            qds[Mapping(qcontig[i][Flipping(j,sz)]*Sgn(j,sz))] = Mapping(qcontig[i][0]);
        }
        qtelo[Mapping(qcontig[i][0])] = qtelo[Mapping(-qcontig[i][sz-1])] = true;
    }
    rep(i,0,r_size) {
        int sz = SZ(rcontig[i]);
        if (sz == 0) continue;
        rep(j,0,2*sz-1) {
            if (j == 0) rpre[Mapping(rcontig[i][0])] = Mapping(-rcontig[i][0]), rnxt[Mapping(-rcontig[i][0])] = Mapping(rcontig[i][0]);
            else rpre[Mapping(rcontig[i][Flipping(j,sz)]*Sgn(j,sz))] = Mapping(rcontig[i][Flipping(j-1,sz)]*Sgn(j-1,sz)), rnxt[Mapping(rcontig[i][Flipping(j-1,sz)]*Sgn(j-1,sz))] = Mapping(rcontig[i][Flipping(j,sz)]*Sgn(j,sz));
            rds[Mapping(rcontig[i][Flipping(j,sz)]*Sgn(j,sz))] = Mapping(rcontig[i][0]);
        }
        rtelo[Mapping(rcontig[i][0])] = rtelo[Mapping(-rcontig[i][sz-1])] = true;
    }
}

Pi findCycle(int x) {
    Pi res = mp(0,0);
    do {
        vis[x] = true;
        if (qtelo[x] || rtelo[x]) res.fi ? (res.se = x) : (res.fi = x);
        x = rnxt[qpre[x]];
    }while (vis[x] == false);
    return res;
}

int findQSet(int x){return qds[x] = (qds[x] == x ? x : findQSet(qds[x]));}

int findRSet(int x){return rds[x] = (rds[x] == x ? x : findQSet(rds[x]));}

void unionQPath(int x, int y) {
    int xpre = qpre[x], ypre = qpre[y];
    qpre[y] = xpre; qpre[x] = ypre; qnxt[xpre] = y; qnxt[ypre] = x;
    qtelo[x] = qtelo[y] = false;
    qds[findQSet(x)] = findQSet(y);
}

void unionRPath(int x, int y) {
    int xpre = rpre[x], ypre = rpre[y];
    rpre[y] = xpre; rpre[x] = ypre; rnxt[xpre] = y; rnxt[ypre] = x;
    rtelo[x] = rtelo[y] = false;
    rds[findRSet(x)] = findRSet(y);
}

void reverseContig(string & s) {
    int sz = SZ(s);
    rep(i,0,sz-1) {
        if (s[i] == 'A') s[i] = 'T';
        else if (s[i] == 'T') s[i] = 'A';
        else if (s[i] == 'C') s[i] = 'G';
        else if (s[i] == 'G') s[i] = 'C';
        else if (s[i] == 'a') s[i] = 't';
        else if (s[i] == 't') s[i] = 'a';
        else if (s[i] == 'c') s[i] = 'g';
        else if (s[i] == 'g') s[i] = 'c';
    }
    reverse(all(s));
}


int main (int argc, char * argv[])
{
// generate .coords file
    checkFile(argc, argv);
    cleanContig();
    string cmd1 = "nucmer " + rcontig_file + " " + qcontig_file + " -p " + temp_folder + "/cssseq";
    string cmd2 = "delta-filter -r -q " + temp_folder + "/cssseq.delta > " + temp_folder + "/cssseq.filter";
    string cmd3 = "show-coords -l -d "  + temp_folder + "/cssseq.filter > " + temp_folder + "/cssseq.coords";
    system(cmd1.c_str()); system(cmd2.c_str()); system(cmd3.c_str());
// represent contig  as a linear permutation
    buildContig();
// make the most rings in the adjacency graph
    initLS();
    vis.resize(2*gene_num+1, false);
    rep(i,1,2*gene_num) {
        if (vis[i]) continue;
        Pi tmp = findCycle(i);
        if (tmp.fi == tmp.se) ;
        else if (qtelo[tmp.fi] && qtelo[tmp.se]) {
            if (findQSet(tmp.fi) != findQSet(tmp.se)) unionQPath(tmp.fi, tmp.se);
        }
        else if (rtelo[tmp.fi] && rtelo[tmp.se]) {
            if (findRSet(tmp.fi) != findRSet(tmp.se)) unionRPath(tmp.fi, tmp.se);
        }
    }
// Output scaffold file
    ifstream input(qcontig_file.c_str());
    ofstream output((output_folder+"/scaffold.fasta").c_str());
    string line;
    while (getline(input, line)) {
        if (line[0] == '>') continue;
        contig.pb(line);
    }
    input.close();
    int scaffold_num = 0, g_num = 0, sum_len = 0;
    Vi temp;
    rep(i,1,2*gene_num) {
        if (qtelo[i] == false) continue;
        output << ">Scaffold_" << ++scaffold_num << endl;
        qtelo[i] = false;
        int x = i;
        do {
            if (gene[abs(AntiMapping(x))-1].r1 == INF) line = contig[gene[abs(AntiMapping(x))-1].id].substr(gene[abs(AntiMapping(x))-1].l1-1);
            else line = contig[gene[abs(AntiMapping(x))-1].id].substr(gene[abs(AntiMapping(x))-1].l1-1, gene[abs(AntiMapping(x))-1].r1-gene[abs(AntiMapping(x))-1].l1+1);
            if (AntiMapping(x) < 0) reverseContig(line);
            output << line;
            temp.pb(abs(AntiMapping(x)));
            sum_len += SZ(line);
            x = qnxt[x];
            ++g_num;
        }while (qtelo[x] == false);
        qtelo[x] = false;
        output << endl;
    }
    rep(i,0,q_size) if (SZ(qcontig[i]) == 0) output << ">Scaffold_" << ++scaffold_num << endl << contig[i] << endl;
    output.close();
    return 0;
}
