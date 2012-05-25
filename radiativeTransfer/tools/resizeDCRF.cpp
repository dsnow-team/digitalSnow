#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cstring>

using namespace std;

//fonction qui prend en entree un fichier brdf produit par pbrt version "photon" et en renvoie un propre a etre trace sous gnuplot


int main(int argc, char *argv[]){

//deltaAngle est la precision que l'on aura sur theta et phi les coordonnees spheriques
int i;
float deltaAngle=15;
string fichier_entree;

if (argc < 2) {cout <<" syntax is < command > filename.txt [--dAngle int]"; return 0;}

for (i=1; i<argc;i++){
	if (!strcmp(argv[i],"--dAngle")) deltaAngle=atof(argv[++i]);
	else if (!strcmp(argv[i],"--help")){ cout <<" syntax is < command > filename.txt [--dAngle int]"; return 0;}
	else fichier_entree=argv[i];
}

deltaAngle=floor(deltaAngle);

if (deltaAngle<1) {
	deltaAngle=1;
	cout <<"la precision sur l'angle est de 1 degre\n";
}
else if (deltaAngle > 90) 
{
	deltaAngle=90;
	cout <<"la precision sur l'angle est de 90 degre\n";
}

string fichier_sortie;
int j;

//on initialise les fichiers entree et sortie
size_t pos=fichier_entree.find(".");
if (pos!=string::npos)
fichier_sortie=fichier_entree.substr(0,pos);
else fichier_sortie=fichier_entree;
fichier_sortie+="ResizeDCRF.txt";

ifstream fichierEntree(fichier_entree.c_str());
ofstream fichierSortie(fichier_sortie.c_str());

//theta et phi sont les coordonnes polaires, gain est un tableau dont la taille depend de deltaAngle 
float theta(0), phi(0), nombrePhoton(0);
int nLignes(floor(90/deltaAngle));
int nColonnes(4*nLignes);
float gain[nLignes][nColonnes];
int phiDegre(0), thetaDegre(0);
int nombrePhotonTotal(0);
char b;
string a,ligne;

for (i=0;i<nLignes;i++)
	for (j=0;j<nColonnes; j++)
		gain[i][j]=0;


fichierEntree.seekg(0,ios::beg);

//on elimine les lignes de commentaires en tete de fichier
while (true)
{	
	fichierEntree >> a ;
	if (isdigit(a.c_str()[0])) break;
	else if (a.c_str()[0]=='#')
		{
		fichierEntree.get(b);
		if (b!='\n') getline(fichierEntree,ligne);
		}
}


//on remplit le tableau gain avec la premiere ligne du fichier
theta=atof(a.c_str());
fichierEntree >>phi;
fichierEntree >> nombrePhoton;
phiDegre=floor(phi);
thetaDegre=floor(theta); 		
if (phiDegre>=90) phiDegre=89;
if (phiDegre<0) phiDegre=0;	
if (thetaDegre <0) thetaDegre=0;		
if (thetaDegre >= 360) thetaDegre=359;
if (thetaDegre < 0) thetaDegre=0;
gain[phiDegre*nLignes/90][thetaDegre*nLignes/90]+=nombrePhoton;
nombrePhotonTotal+=nombrePhoton;

//on remplit gain avec le reste du fichier
while (!fichierEntree.eof())
{
	fichierEntree >> theta;
	fichierEntree >> phi;
	fichierEntree >> nombrePhoton;

	phiDegre=floor(phi);
	thetaDegre=floor(theta); 			
	if (phiDegre>=90) phiDegre=89;
	if (phiDegre<0) phiDegre=0;	
	if (thetaDegre <0) thetaDegre=0;		
	if (thetaDegre >= 360) thetaDegre=359;
	if (thetaDegre < 0) thetaDegre=0;
	gain[phiDegre*nLignes/90][thetaDegre*nLignes/90]+=nombrePhoton;
	nombrePhotonTotal+=nombrePhoton;
}



//on normalise le gain pour avoir la DCRF
for (i=0;i<nLignes;i++)
	for (j=0;j<nColonnes;j++)
		gain[i][j]=gain[i][j]*180/(M_PI*deltaAngle*cos(i*M_PI*deltaAngle/180)*fabs((cos((i+1)*M_PI*deltaAngle/180)-cos(i*M_PI*deltaAngle/180))));
		


//on renvoie dans le fichier de sortie les valeurs
fichierSortie <<"# this file provides the anisotropic reflectance factor R, a variant of the BRDF :\n# R = BRDF * Pi / Albedo \n# pour tracer avec gnuplot :\n#splot \"file.txt\" using 1:2:3:4 with pm3d\n"; 
	
for (j=0;j<360; j++){
	for (i=0;i<90;i++)
		{fichierSortie << cos(j*M_PI/180)*cos((90-i)*M_PI/180) << " " << sin(j*M_PI/180)*cos((90-i)*M_PI/180) << " " << sin((90-i)*M_PI/180) << " " << gain[i*nLignes/90][j*nLignes/90]*M_PI/nombrePhotonTotal<< endl; 					
		}

	fichierSortie << endl;		
}

for (i=0;i<90;i++)
	fichierSortie << cos(360*M_PI/180)*cos((90-i)*M_PI/180) <<" " << sin(360*M_PI/180)*cos((90-i)*M_PI/180) << " " << sin((90-i)*M_PI/180) << " " << gain[i*nLignes/90][0]*M_PI/nombrePhotonTotal<< endl; 		

return 0;

}
