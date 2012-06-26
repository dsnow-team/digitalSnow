#include <fstream>
#include <string>
#include <cstdlib>
#include <iostream>
#include <cctype>
#include <cstring>
#include <algorithm>


using namespace std;
//variables globales qui cernent la bounding box
double minX(0), maxX(0), minY(0), maxY(0), minZ(0), maxZ(0);

void ecritFichierGeometrie(string fichierNoff, string fichierGeomPbrt);

void ecritFichierPbrt(string fichierPbrt, string fichierGeomPbrt, string fichierEXR);

void ecritFichierPhoton(string fichierPhoton, string fichierGeomPbrt);



int main(int argc, char *argv[])
{

  cout << "Warning !! the file provided must be noff without blank lines and with commentary only at the top !\nThe normals provided MUST be inwarded"<< endl; 
  if (argc < 2) {
    cout<< "bad syntax. The syntax is : <command> -i input.off -o output";
    exit(2);
  }

  string fichierNoff, fichier_sortie;
  bool entre(false), sortie(false);


  for (int i=1; i<argc;i++){
    if (!strcmp(argv[i],"--help") || !strcmp(argv[i],"--help")){cout << "syntax : <command> -i input.noff -o output\n"; return 0;}
    else if (!strcmp(argv[i],"--input") || !strcmp(argv[i],"-i")) {fichierNoff=argv[++i]; entre=true;}
    else if (!strcmp(argv[i],"--output") || !strcmp(argv[i],"-o")) {fichier_sortie=argv[++i]; sortie=true;}
  }

  if (!entre || !sortie) 
    {
      cout << "syntax : <command> -i input.noff -o output\n"; 
      exit(1);
    }
  //on prend en entrée un fichier noff et on sort 2 fichier : un de geometrie et le corps du fichier .pbrt


  string fichierGeomPbrt,fichierPbrt,fichierPhoton, fichierEXR;

  fichierGeomPbrt=fichier_sortie;
  fichierGeomPbrt+="Geometry.pbrt";
  fichierPhoton=fichier_sortie;
  fichierPhoton+="Photon.pbrt";
  fichierPbrt=fichier_sortie;
  fichierPbrt+="Image.pbrt";
  fichierEXR=fichier_sortie;
  fichierEXR+=".exr";


  ecritFichierGeometrie(fichierNoff, fichierGeomPbrt);

  ecritFichierPbrt(fichierPbrt,fichierGeomPbrt,fichierEXR);

  ecritFichierPhoton(fichierPhoton, fichierGeomPbrt);

  cout <<"the length of the image file is "<<maxX-minX <<" * " << maxY -minY << " * "<< maxZ-minZ <<endl;

  return 0;
}


//la fonction qui ecrit le fichier de geometrie

void ecritFichierGeometrie(string fichierNoff, string fichierGeomPbrt)
{


  ifstream fichierEntree(fichierNoff.c_str());
  ofstream fichierSortieGeom(fichierGeomPbrt.c_str());



  //on a 2 fichiers temporaires qui contiennent les points et les vecteurs normaux
  ofstream fichierPoint("pointTemp.txt");
  ofstream fichierVecteur("vecteurTemp.txt");

  //initialisation des variables : nombrePoints et nombreFaces sont explicites
  //nombreVertex sert pour compter le nombre de sommets d'une face
  //min et max servent à calculer la bounding box
  string ligne;
  int i(0),j(0),nombrePoints(0), nombreFaces(0), nombreVertex(0);
  size_t pos1;
  string a;
  char b;
  int indice(0), indice1(0);
  float point[3], vecteur[3];

  //on se place à l'entrée du fichier et on saute les commentaires

  fichierEntree.seekg(0,ios::beg);
  while (true)
    {	
      fichierEntree >> a;
      if (isdigit(a.c_str()[0])) break;
      else if (a.c_str()[0]=='#')
	{
	  fichierEntree.get(b);
	  if (b!='\n') getline(fichierEntree,ligne);
	}
    }


  // on initialise le nombre de points et le nombre de faces
  nombrePoints=atoi(a.c_str());
  fichierEntree >> nombreFaces;
  getline(fichierEntree,ligne);

 
  //initialisation des min et max
  fichierEntree >> point[0];
  minX=point[0];
  maxX=point[0];	


  fichierEntree >> point[1];
  minY=point[1];
  maxY=point[1];	


  fichierEntree >> point[2];
  minZ=point[2];
  maxZ=point[2];	


  fichierEntree >> vecteur[0];
  fichierEntree >> vecteur[1];
  fichierEntree >> vecteur[2];

  //on commence à remplir les fichiers point et vecteur

  fichierPoint << point[0] << " ";
  fichierPoint << point[1]<< " ";
  fichierPoint << point[2] << "\n";


  fichierVecteur << -vecteur[0] <<" " ;

  fichierVecteur << -vecteur[1] << " ";

  fichierVecteur << -vecteur[2] <<"\n";


  // on remplit les fichiers points et vecteur jusqu'à ce qu'il n'y ait plus de points
  for (i=1;i<nombrePoints;i++)
    {
	
      fichierEntree >> point[0];
      if (point[0]<minX) minX=point[0];	
      if (point[0]>maxX) maxX=point[0];	
      fichierPoint << point[0] << " "; 

      fichierEntree >> point[1];
      if (point[1]<minY) minY=point[1];	
      if (point[1]>maxY) maxY=point[1];		
      fichierPoint << point[1]<< " "; 

      fichierEntree >> point[2];
      if (point[2]<minZ) minZ=point[2];	
      if (point[2]>maxZ) maxZ=point[2];	
      fichierPoint << point[2] << "\n"; 
	
      fichierEntree >> vecteur[0];
      fichierEntree >> vecteur[1];
      fichierEntree >> vecteur[2];
      fichierVecteur << -vecteur[0] << " ";
      fichierVecteur << -vecteur[1] << " ";
      fichierVecteur << -vecteur[2] << "\n";

    }

  fichierPoint.close();
  fichierVecteur.close();


  //on remplit le fichier de geométrie avec les points et vecteurs

  fichierSortieGeom << "Shape \"trianglemesh\" \"point P\" [ \n";

  ifstream fichierpoint("pointTemp.txt");
  while (getline(fichierpoint, ligne)){
    fichierSortieGeom << ligne;
    fichierSortieGeom << "\n";
  }



  fichierSortieGeom << "] \"normal N\" [\n";
  ifstream fichiervect("vecteurTemp.txt");
  while (getline(fichiervect, ligne)){
    fichierSortieGeom << ligne;
    fichierSortieGeom << "\n";
  }


  //on supprime les fichiers temporaires
  remove("pointTemp.txt");
  remove("vecteurTemp.txt");

  //on écrit les indices des faces
  fichierSortieGeom << "] \"integer indices\" ["; 

  for (i=0;i<nombreFaces ; i++)
    {
      fichierEntree >> nombreVertex;
      fichierEntree >> indice;
      fichierEntree >> indice1;
	
      for (j=0; j<nombreVertex -2;j++){
	fichierSortieGeom << indice << " " <<indice1 << " ";	
	fichierEntree >> indice1;
	fichierSortieGeom << indice1 << "\n";
      }
    }
  fichierSortieGeom << "]";

  cout << "geometry file has been released"<<endl; 

}






//la fonction qui sort le fichier pbrt lisible par le logiciel
void ecritFichierPbrt(string fichierPbrt, string fichierGeomPbrt, string fichierEXR){

  ofstream fichierSortiePbrt(fichierPbrt.c_str());
double Maximum(0);

Maximum = max(max(maxX-minX, maxY-minY), maxZ-minZ);
if (Maximum==0){cout << "geometry file is empty !!!" <<endl; Maximum=1;} 


  // et on écrit dans le fichier qu'il faudra lancer sous pbrt
  fichierSortiePbrt <<"## to have a nicest image level up the \"samplesperpixel\" of metropolis\n## to be more rapid, make this number down (but loose quality of image)\n \n \n";

  //declaration des attributs generaux
  fichierSortiePbrt << "Scale -1.000000 1.000000 1.000000 \n \nTranslate -278.000000 -273.000000 500.000000\n \nRenderer \"metropolis\" \"integer samplesperpixel\" [128]\n \nCamera \"perspective\" \"float fov\" [55.000000]\n \nFilm \"image\" \"integer xresolution\" [1000] \"integer yresolution\" [750]\n    \"string filename\" \""<< fichierEXR  <<"\"\n \nPixelFilter \"box\" \n \nWorldBegin\n \n AttributeBegin\nTranslate 340.000000 278.000000 -50\nLightSource \"point\" \"point from\" [0.000000 200.000000 -50.000000] \"color I\" [412300 341100 298600]\nAttributeEnd\n \n";

  //declaration du fond
  fichierSortiePbrt << "#le fond \nAttributeBegin\nMaterial \"matte\" \"color Kd\" [.8 .8 .8]\nShape \"trianglemesh\"  \"integer indices\" [0 2 1 0 3 2] \"point P\" [800.000000 0.000000 0.000000 -250.000000 0.000000 0.000000 -250.000000 0.000000 1000.000000 800.000000 0.000000 1000.000000]\nAttributeEnd\n \n";

  //on inclut l'echantillon
  fichierSortiePbrt << "#include of the sample\nAttributeBegin\nTranslate 180 180 386 \nRotate -30 1 0 0\nRotate 30 0 1 0\nScale "<< (double)256/Maximum << " " << (double)256/Maximum<<" " << (double)256/Maximum<<"\n Rotate -90 1 0 0 \nTranslate "<< -minX <<" " << -minY << " " <<-minZ <<"\nMaterial \"matte\" \"color Kd\" [1 1 1]\nInclude \"" << fichierGeomPbrt<<"\"\nAttributeEnd\n \nWorldEnd";

  cout << "pbrt file generated"<<endl;
}




//la fonction qui sort le fichier pour le lanceur de photon
void ecritFichierPhoton(string fichierPhoton, string fichierGeomPbrt){
  double facteur=(maxZ-minZ)/256;

  ofstream fichierSortiePhoton(fichierPhoton.c_str());

  //on definit la photonmap
  fichierSortiePhoton << "## causticphotons = number of launched photons;\n## maxdepth= max number of intersections for one photon before stopping\n\nSurfaceIntegrator \"photonmap\" \"integer indirectphotons\" [0] \"integer causticphotons\" [20000]\n\"integer maxspeculardepth\" [100000] \"integer maxphotondepth\" [100000]\n\n";

  //definition de la source de lumiere
  fichierSortiePhoton << "#ligth source : to change the direction of the source change point from and point to, only the direction is important \nWorldBegin\n \nAttributeBegin\nLightSource \"distant\" \"point from\" [0 0 50] \"point to\" [0 0 0]\nAttributeEnd\n\n";

  //pour dupliquer l'echantillon et prendre en compte la profondeur
  fichierSortiePhoton << "#to duplicate the sample and take care of the depth, we create a cube wich will surround the sample\n\nAttributeBegin\nMaterial \"glass\"\nShape \"trianglemesh\" \"integer indices\" [0 1 2]\n\"point P\" [0 0 256.001 " << (maxX-minX)/facteur << " 0 256.001 0 " << (maxY-minY)/facteur << " 256.001]\n\"normal N\" [0 0 1 0 0 1 0 0 1]\n\nShape \"trianglemesh\" \"integer indices\" [0 1 2]\n\"point P\" [ "<< (maxX-minX)/facteur<< " 0 256.001 0 "<<(maxY-minY)/facteur<< " 256.001 "<<(maxX-minX)/facteur<< " "<<(maxY-minY)/facteur<< " 256.001]\n\"normal N\" [0 0 1 0 0 1 0 0 1]\n\nShape \"trianglemesh\" \"integer indices\" [0 1 2]\n\"point P\" [0 0 0 " << (maxX-minX)/facteur <<" 0 0 0 " << (maxY-minY)/facteur << " 0]\n\"normal N\" [0 0 -1 0 0 -1 0 0 -1]\n\nShape \"trianglemesh\" \"integer indices\" [0 1 2]\n\"point P\" [" << (maxX-minX)/facteur <<" 0 0 0 " << (maxY-minY)/facteur <<" 0 " << (maxX-minX)/facteur <<" " << (maxY-minY)/facteur <<" 0]\n\"normal N\" [0 0 -1 0 0 -1 0 0 -1]\n\nShape \"trianglemesh\" \"integer indices\" [0 2 1]\n\"point P\" [0 0 0 " << (maxX-minX)/facteur <<" 0 0 0 0 256]\n\"normal N\" [0 -1 0 0 -1 0 0 -1 0]\n\nShape \"trianglemesh\" \"integer indices\" [0 2 1]\n\"point P\" [" << (maxX-minX)/facteur <<" 0 0 0 0 256 " << (maxX-minX)/facteur <<" 0 256]\n\"normal N\" [0 -1 0 0 -1 0 0 -1 0]\n\nShape \"trianglemesh\" \"integer indices\" [0 1 2]\n\"point P\" [0 " << (maxY-minY)/facteur <<" 0 " << (maxX-minX)/facteur <<" " << (maxY-minY)/facteur <<" 0 0 " << (maxY-minY)/facteur <<" 256]\n\"normal N\" [0 1 0 0 1 0 0 1 0]\n\nShape \"trianglemesh\" \"integer indices\" [0 2 1]\n\"point P\" [" << (maxX-minX)/facteur <<" " << (maxY-minY)/facteur <<" 0 0 " << (maxY-minY)/facteur <<" 256 " << (maxX-minX)/facteur <<" " << (maxY-minY)/facteur <<" 256]\n\"normal N\" [0 1 0 0 1 0 0 1 0]\n\nShape \"trianglemesh\" \"integer indices\" [0 1 2]\n\"point P\" [" << (maxX-minX)/facteur <<" 0 256 " << (maxX-minX)/facteur <<" 0 0 " << (maxX-minX)/facteur <<" " << (maxY-minY)/facteur <<" 0]\n\"normal N\" [1 0 0 1 0 0 1 0 0]\n\nShape \"trianglemesh\" \"integer indices\" [0 1 2]\n\"point P\" [" << (maxX-minX)/facteur <<" " << (maxY-minY)/facteur <<" 0 " << (maxX-minX)/facteur <<" 0 256 " << (maxX-minX)/facteur <<" " << (maxY-minY)/facteur <<" 256]\n\"normal N\" [1 0 0 1 0 0 1 0 0]\n\nShape \"trianglemesh\" \"integer indices\" [0 1 2]\n\"point P\" [0 0 0 0 0 256 0 " << (maxY-minY)/facteur <<" 0]\n\"normal N\" [-1 0 0 -1 0 0 -1 0 0]\n\nShape \"trianglemesh\" \"integer indices\" [0 1 2]\n\"point P\" [0 0 256 0 " << (maxY-minY)/facteur <<" 0 0 " << (maxY-minY)/facteur <<" 256]\n\"normal N\" [-1 0 0 -1 0 0 -1 0 0]\n\nAttributeEnd\n";

  //on inclut le fichier de geometrie
  fichierSortiePhoton <<"\n#path to the geometry file \nAttributeBegin\nMaterial \"glass\"\nScale "<<(double)256/(maxZ-minZ)<< " " << (double)256/(maxZ - minZ)<<" " << (double)256/(maxZ-minZ)<<"\nTranslate " <<-minX << " " << -minY <<" " << -minZ <<"\nInclude \""<< fichierGeomPbrt <<"\"\nAttributeEnd\n\nWorldEnd";

  cout <<"photon pbrt file generated" <<endl;


}
