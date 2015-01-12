#include <fstream>
#include <new>

//display 3D
// static 
  #ifdef WITH_CAIRO
#include "DGtal/io/boards/Board3DTo2D.h"
  #endif

#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/writers/VolWriter.h"

#include "LocalMCM.h"

#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/kernel/BasicPointPredicates.h"

// Cell embedders
#include <DGtal/topology/CanonicCellEmbedder.h>
#include <DGtal/topology/CanonicSCellEmbedder.h>

// frontier
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/topology/helpers/FrontierPredicate.h"
#include "DGtal/topology/LightExplicitDigitalSurface.h"

// interactive
#include <QtGui/qapplication.h>
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/CDrawableWithDisplay3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"

using namespace DGtal::functors;

template< typename TViewer, typename TImage >
bool displayPartition(TViewer& viewer, const TImage& img)
{
  typedef typename TImage::Value Label; 

  //KhalimskySpace
  Domain d = img.domain();
  Point aLowerBound = d.lowerBound(); 
  Point aUpperBound = d.upperBound(); 
  KSpace aKSpace;
  aKSpace.init(aLowerBound, aUpperBound, true);

  // Update embedders for the viewer
  viewer.setKSpaceEmbedder(  new CanonicCellEmbedder<KSpace>(aKSpace) );
  viewer.setSKSpaceEmbedder( new CanonicSCellEmbedder<KSpace>(aKSpace) );

  //container
  std::set<Cell> aSet; 

  //mark bels
  for (DGtal::Dimension k = 0; k < KSpace::dimension; ++k )
    {
      Cell dir_low_uid = aKSpace.uSpel( aLowerBound );
      Cell dir_up_uid = aKSpace.uGetDecr( aKSpace.uSpel( aUpperBound ), k);
      Cell p = dir_low_uid;
      do 
        {
          Label here = img( aKSpace.uCoords(p) );
          Label next = img( aKSpace.uCoords(aKSpace.uGetIncr( p, k )) );
          if ( here != next ) 
            { // add new bel to the set.
              aSet.insert( aKSpace.uIncident( p, k, true ));
            }
        }
      while ( aKSpace.uNext( p, dir_low_uid, dir_up_uid ) );
    }

  GradientColorMap<long> colorMap( 0, 510 );
  colorMap.addColor(Color::Yellow);
  colorMap.addColor(Color::Blue);
  colorMap.addColor(Color::Red);
  colorMap.addColor(Color::Green);

  /// retrieve frontiers
  unsigned int counter = 0; 
  while(!aSet.empty()){
 
    SCell sbel = aKSpace.signs( *(aSet.begin()), true ); 
    //incident points
    SCellToIncidentPoints<KSpace> func( aKSpace ); 
    typename SCellToIncidentPoints<KSpace>::Output points = func( sbel ); 
    Label iLabel( img( points.first ) ); 
    Label oLabel( img( points.second ) ); 

    /// frontier from sbel
    typedef FrontierPredicate<KSpace, TImage> SurfelPredicate;
    /// !!!!!! be careful oLabel and iLabel are swaped because func is wrong
    ///     => Seems to have been corrected since ...
    SurfelPredicate surfelPred( aKSpace, img, iLabel, oLabel ); 
    typedef LightExplicitDigitalSurface<KSpace, SurfelPredicate> Frontier;
    Frontier frontier( aKSpace, 
		       surfelPred, 
		       SurfelAdjacency<KSpace::dimension>( true ), 
		       sbel ); 

    // tracking (and removing bels belonging to this frontier)
    // and display
    counter++; 
    typedef typename Frontier::SurfelConstIterator SurfelIterator;
    for ( SurfelIterator it = frontier.begin(), 
	    itEnd = frontier.end();
	  it != itEnd; ++it )
      {
	viewer << DGtal::CustomColors3D( colorMap( iLabel+oLabel ), 
					 colorMap( iLabel+oLabel ) );
	viewer << aKSpace.unsigns( *it );
	aSet.erase( aKSpace.unsigns( *it ) );
      }
  }
  trace.info() << counter << " frontier(s) displayed" << std::endl; 

  // What should it return ?
  return true;
}

template< typename TImage >
bool writePartition(const TImage& img, string filename, string format)
{

  if (format.compare("pngc")==0)
  {
  #ifdef WITH_CAIRO
    Board3DTo2D<> viewer;
    
    displayPartition( viewer, img ); 
    // Domain d = img.domain(); 
    // Domain::ConstIterator cIt = d.begin(); 
    // Domain::ConstIterator cItEnd = d.end(); 
    // for ( ; cIt != cItEnd; ++cIt)
    // { 
    //   if (img(*cIt) <= threshold) 
    // 	      viewer << *cIt; 
    // }

  //reading camera configuration for the 3d to 2d projection
  std::ifstream file(".camera", ios::in); 
  if(file) 
  {       

    std::vector<std::vector<double> > p; //position/direciton
    double znear, zfar; 

    std::string line;
    getline(file, line); //skip the first line
    //for the three following lines
    unsigned int nbLines = 1; 
    while( std::getline(file, line) && (nbLines <= 3) ) 
    {
      p.push_back( std::vector<double>(3) ); 

      std::istringstream isline( line );
      std::string word;
      std::getline( isline, word, ' ' ); //skip the first word
      unsigned int k = 0; 
      while ( std::getline( isline, word, ' ' ) && (k <= 3) )
      {
          std::istringstream isword( word.substr (0,word.size()-1) );
          isword >> p[nbLines-1][k];
          ++k; 
      }
      ++nbLines; 
    }
    //for the last line
    if ( std::getline(file, line) ) 
    { 
      std::istringstream isline( line );
      std::string word;
      std::getline( isline, word, ' ' ); //skip the first word
      std::getline( isline, word, ' ' );
      std::istringstream is1( word );
      is1 >> znear;
      std::getline( isline, word, ' ' ); //skip the third word
      std::getline( isline, word, ' ' ); //skip the fourth word
      std::getline( isline, word, ' ' );
      std::istringstream is2( word );
      is1 >> zfar;
    }
    file.close();

    //setting camera configuration
    viewer << CameraPosition(p[0][0], p[0][1], p[0][2])
	   << CameraDirection(p[1][0], p[1][1], p[1][2]) 
     << CameraUpVector(p[2][0], p[2][1], p[2][2])
     << CameraZNearFar(znear, zfar);

  }
  else  
  {
    trace.emphase() << "Failed to read '.camera'. Default camera configuration" << std::endl;

    //default config
    typename TImage::Vector v = img.extent(); 
    viewer << CameraPosition(v[0]/2, v[1]/2, 2*v[2])
	   << CameraDirection(0, 0, -1) 
     << CameraUpVector(0, 1, 0)
     << CameraZNearFar(v[2]/2, 3*v[2]);
  }

    int size = img.extent()[0]; 
    std::stringstream s; 
    s << filename << ".png";
    viewer.saveCairo(s.str().c_str(),Board3DTo2D<>::CairoPNG,3*size/2,3*size/2 ); 
    return true; 
  #else
    trace.emphase() << "Failed to use Cairo 3d to 2d (not installed)" << std::endl;
    return false; 
  #endif

  } else if (format.compare("png")==0)
    {
    trace.emphase() << "snapshot with QGLViewer" << std::endl;

    int argc = 1; 
    string sargv1 = "QGLViewer"; 
    char* argv1 = const_cast<char*>( sargv1.c_str() ); 
    char* argv[1]; 
    argv[0] = argv1; 
    QApplication application(argc,argv);
    Viewer3D<> viewer;
    viewer.show();

    //display
    displayPartition(viewer, img); 
    viewer << Viewer3D<>::updateDisplay;

    if (QGLViewer::QGLViewerIndex(&viewer) > 0)
      {//rename state file
	string oldf = ".qglviewer.xml";
	std::stringstream news; 
	news << ".qglviewer" << (QGLViewer::QGLViewerIndex(&viewer)) << ".xml";
	string newf = news.str();  
	if (rename (oldf.c_str(), newf.c_str()) == -1) 
	  trace.info() << "renaming " << oldf << " into " 
		       << newf << " failed " << std::endl; 
      }

    if (!viewer.restoreStateFromFile())
      {
	string s = viewer.stateFileName().toStdString(); 
	trace.emphase() << " file " << s 
		      << " not found " 
			<< std::endl;
      }
    viewer.updateGL(); 

    viewer.setSnapshotFileName(filename.c_str());  
    viewer.setSnapshotFormat("PNG");  
    viewer.saveSnapshot(true, true); 

    {//rename snapshot
    std::stringstream olds;
    olds << viewer.snapshotFileName().toStdString()
	 << "-" << setfill('0') << std::setw(4) 
	 << (viewer.snapshotCounter()-1) << ".png"; 
    string oldf = olds.str(); 
    std::stringstream news; 
    news << filename << ".png";
    string newf = news.str();  
    if (rename (oldf.c_str(), newf.c_str()) == -1) 
      trace.info() << "renaming " << oldf << " into " 
		   << newf << " failed " << std::endl; 
    }
    
    {//rename state file
    string oldf = viewer.stateFileName().toStdString();
    std::stringstream news; 
    news << ".qglviewer" << (QGLViewer::QGLViewerIndex(&viewer)+1) << ".xml";
    string newf = news.str();  
    if (rename (oldf.c_str(), newf.c_str()) == -1) 
      trace.info() << "renaming " << oldf << " into " 
		   << newf << " failed " << std::endl; 
    }

    viewer.setStateFileName(QString::null);  
    application.exit();

    // Guess ...
    return true;

    }
  else if (format.compare("vol")==0)
  {

    //write it into a vol file
    std::stringstream s; 
    s << filename << ".vol";
    typedef Cast<unsigned char> Fonctor; 
    VolWriter<TImage, Fonctor>::exportVol( s.str(), img, Fonctor() );

    return true; 

 } else return false; 
  
}



// template< typename TImage >
// bool displayImage(int argc, char** argv, const TImage& img, const double& threshold = 0)
// {

//   //KhalimskySpace
//   Domain d = img.domain(); 
//   KSpace K;
//   K.init(d.lowerBound(), d.upperBound(), true);
//   //adjacency  
//   SurfelAdjacency<3> SAdj( true );
//   std::vector<std::vector<SCell> > vectConnectedSCell;
//   //predicate
//   typedef SimpleThresholdForegroundPredicate<TImage> PointPredicate; 
//   PointPredicate predicate(img,threshold);
//   //tracking 
//   Surfaces<KSpace>::extractAllConnectedSCell(vectConnectedSCell,K, SAdj, predicate, true);

//   #ifdef WITH_VISU3D_QGLVIEWER
//   QApplication application(argc,argv);
//   Viewer3D<> viewer;
//   viewer.show();

//   for(unsigned int i=0; i< vectConnectedSCell.size();i++){
//     for(unsigned int j=0; j< vectConnectedSCell.at(i).size();j++){
//       viewer << vectConnectedSCell.at(i).at(j);
//     }    
//   }

//   viewer << Viewer3D<>::updateDisplay;

//   return application.exec();
// #else
//   return false; 
// #endif
// }



template< typename TLabelImage, typename TDistanceImage, typename TExternImage >
bool displayImageWithInfo(int argc, char** argv, const TLabelImage& limg, 
		   TDistanceImage& img, 
		  const TExternImage& ext1, const TExternImage& ext2, 
		  const short int& threshold = 0)
{

  //KhalimskySpace
  Domain d = img.domain(); 
  KSpace K;
  K.init(d.lowerBound(), d.upperBound(), true);
  //adjacency  
  SurfelAdjacency<Space::dimension> SAdj( true );
  std::vector<std::vector<SCell> > vectConnectedSCell;
  //predicate
  typedef Thresholder<typename TLabelImage::Value,true,true> Binarizer; 
  Binarizer b(threshold); 
  PointFunctorPredicate<TLabelImage,Binarizer> predicate(limg, b);
  //tracking 
  Surfaces<KSpace>::extractAllConnectedSCell(vectConnectedSCell,K, SAdj, predicate, true);
  //NB. tracking is done on the label image because the distance image
  //may be not defined everywhere with a local evolution method

  //local MCM operator 
  typedef LocalMCM<TDistanceImage, TExternImage> DiffOperator; 
  DiffOperator op(img, ext1, ext2); 

  #ifdef WITH_VISU3D_QGLVIEWER
  QApplication application(argc,argv);
 
  Viewer3D<> viewer(K);
  viewer.show();

  //good for Al
  GradientColorMap<double> colorMap4Pos( 0.0, 2.5 );
  colorMap4Pos.addColor( Color( 255, 255, 255 ) );
  colorMap4Pos.addColor( Color( 0, 0, 255 ) );
  GradientColorMap<double> colorMap4Neg( -2.5, 0.0 );
  colorMap4Neg.addColor( Color( 0, 255, 0 ) );
  colorMap4Neg.addColor( Color( 255, 255, 255 ) );

  for(unsigned int i=0; i< vectConnectedSCell.size();i++){
    for(unsigned int j=0; j< vectConnectedSCell.at(i).size();j++){

	SCell s = vectConnectedSCell.at(i).at(j); 
	Point p = K.sCoords( K.sDirectIncident( s, *K.sOrthDirs( s ) ) ); 

	//curvature
	typename DiffOperator::Curvature curvature = op.getCurvature( p );
	if (curvature >= 0)
	    viewer << DGtal::CustomColors3D( colorMap4Pos( curvature ), 
					     colorMap4Pos( curvature ) );
	else 
	    viewer << DGtal::CustomColors3D( colorMap4Neg( curvature ), 
					     colorMap4Neg( curvature ) );
	viewer << s; 


	//naive normal
	//Vector normal = K.sKCoords( s ) - K.sKCoords( K.sDirectIncident( s, *K.sOrthDirs( s ) ) ); 
	//normal
	typename DiffOperator::Normal normal = op.getNormal( p );
	Space::RealPoint center( K.sKCoords(s) );
	center /= 2;
	center -= Space::RealVector(0.5, 0.5, 0.5); 
	double normalNorm = std::sqrt( normal[0]*normal[0] 
				       + normal[1]*normal[1] 
				       + normal[2]*normal[2] );
	normal /= normalNorm; 
    
    viewer.setLineColor( DGtal::Color(255,0,0) );
    viewer.addLine(
        Space::RealPoint(center[0], center[1], center[2]),
        Space::RealPoint(center[0]+normal[0], center[1]+normal[1], center[2]+normal[2]),
        1.0
    );

      }
  }

  viewer << Viewer3D<>::updateDisplay;

  return application.exec();
#else
  return false; 
#endif
}


template< typename TImage >
bool displayPartition(int argc, char** argv, const TImage& img)
{

  #ifdef WITH_VISU3D_QGLVIEWER
  QApplication application(argc,argv);
  Viewer3D<> viewer;
  viewer.show();

  displayPartition(viewer, img); 

  viewer.setSnapshotFormat("PNG");  
  viewer << Viewer3D<>::updateDisplay;

  return application.exec();
#else
  return false; 
#endif
}
