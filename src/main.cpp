#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <csignal>

#ifdef SFML
#include <SFML/System.hpp>
#include <SFML/Window.hpp>
#include <SFML/OpenGL.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/Config.hpp>
#endif

#include "World.h"
#include "Vector3D.h"
#include "NormalDistribution.h"

TYPE_FLOAT ** outputArray; ///< buffer for file output
int numLines; ///< number of lines to write at once
int numEntries; ///< entries per line

#ifdef GL
void glInit();
#endif

void processSIGTERM(int sig)
{
  std::cout<< "SIGTERM RECEIVED!" << std::endl;
  WorldSettings::terminate = true;
}

void processSIGCheckpoint(int sig)
{
  std::cout<< "SIGCheckpoint RECEIVED!" << std::endl;
  WorldSettings::saveCheckpoint = true;
}

void writeOutput(std::ofstream * _file, int currentLine)
{
  std::stringstream outputBlock;
  for (int iLine = 0 ; iLine < currentLine ; iLine++){
    for (int iEntry = 0 ; iEntry < numEntries ; iEntry++){
      outputBlock<<outputArray[iLine][iEntry]<<",";
    }
    outputBlock<<std::endl;
  }
  (*_file) << outputBlock.str();
}

/// \todo check that references are correct
void writeOutputBinary(std::ofstream * _file, int currentLine){
  (*_file).write((char *) &(outputArray[0][0]), currentLine*numEntries*sizeof(TYPE_FLOAT));
}
		      

int main ( int argc, char** argv )
{ 
  std::cout<<time(NULL)<<std::endl;
  std::cout<<sizeof(Vector3D)<<std::endl;

  signal(SIGTERM, processSIGTERM);

  signal(SIGINT, processSIGCheckpoint);

  // file write array
  int currentLine = 0;
  numEntries = 13;
  numLines = 20;
  outputArray = new TYPE_FLOAT*[numLines];
  for(int i = 0; i < numLines; ++i)
    {
      outputArray[i] = new TYPE_FLOAT[numEntries];
    }
  
  int isFinished = 0;

  // check if a different configuration file has been specified
  if (argc == 1)
    WorldSettings::initialize();
  else
    for (int i = 1; i < argc; i+=2) /// \todo why inc by 2?
      WorldSettings::initialize(argv[i]);

  TYPE_FLOAT realDT = WorldSettings::dt;
  World w;
    
  TYPE_FLOAT lastTime = w.getTime();

  if (WorldSettings::outSetting == WorldSettings::FILE)
    {
      /* open the file for output */
      std::stringstream filenamestream;
      filenamestream << WorldSettings::outFile << ".txt";
      std::string filename = filenamestream.str();
      std::ofstream file;
      file.open(filename.c_str(), std::fstream::out | std::fstream::app );
      //file.open(filename.c_str(), std::ios::out | std::ios::app | std::ios::binary);
      file<<"#seed: "<<WorldSettings::seed<<std::endl;
      
      std::cout<<"Writing to: "<<filename<<std::endl;
      
      TYPE_FLOAT lastTime = w.getTime();

      while (w.getStep() < WorldSettings::steps && isFinished != 1)
        {
	  //w.resetDisplacementSquared();
	  realDT = w.simulate();

	  // make a snapshot of current configuration
	  if(WorldSettings::stepsBetweenSnapshot > 0 &&
	     w.getStep() % WorldSettings::stepsBetweenSnapshot ==0){
	    w.savePositions(); 
	    w.saveStateBinaryFast(WorldSettings::snapshotFile);
	  }
	  
	  if(w.getStep() % WorldSettings::stepsBetweenFileOutput == 0)
	    {
	      if (WorldSettings::stepsBetweenSeparation > 0)
		w.updateSeparationIJ(realDT);
	      
	      Vector3D com = w.centreOfMass() - WorldSettings::offsetVector;
	      Vector3D re = w.endToEndVector();
	      Vector3D eigen_rg = w.asphericity();
	 
	      // write the output array to file and reset the counter
	      if ( currentLine == numLines) {
		writeOutput(&file, currentLine);
		currentLine = 0;
	      }
	      
	      outputArray[currentLine][0] = w.getTime();
	      //hydrodynamic radius of the first partile
	      outputArray[currentLine][1] = w.displacementSquared(0)/realDT;
	      outputArray[currentLine][2] = w.rkInverse();
	      outputArray[currentLine][3] = com.x;
	      outputArray[currentLine][4] = com.y;
	      outputArray[currentLine][5] = com.z;
	      outputArray[currentLine][6] = re.x;
	      outputArray[currentLine][7] = re.y;
	      outputArray[currentLine][8] = re.z;
	      outputArray[currentLine][9] = eigen_rg.x;
	      outputArray[currentLine][10] = eigen_rg.y;
	      outputArray[currentLine][11] = eigen_rg.z;
	      outputArray[currentLine][12] = w.averageBondLength();

	      currentLine++;  
		  
	    }

	  if(WorldSettings::stepsBetweenSeparation > 0 &&
	     w.getStep() % WorldSettings::stepsBetweenSeparation ==0){
	    w.saveSeparationIJ(WorldSettings::separationFile);
	  }
	  
	  lastTime = w.getTime();
	  
	  if(WorldSettings::terminate)
	    {
	      //w.saveContactMap(WorldSettings::contactMapFile);
	      w.saveSeparationIJ(WorldSettings::separationFile);
	      w.savePositions();
	      w.saveStateBinary(WorldSettings::checkpointFile);
	      isFinished = 1;
	    }
          if(WorldSettings::saveCheckpoint)
	    {
	      //w.saveContactMap(WorldSettings::contactMapFile);
	      w.savePositions();
	      w.saveStateBinary(WorldSettings::checkpointFile);
	      WorldSettings::saveCheckpoint = false;
	    }
	  
        }
      // write the incomplete output array to file
      writeOutput(&file, currentLine);
      file.close();
    } 
  else if (WorldSettings::outSetting == WorldSettings::SCREEN)
    {
      TYPE_FLOAT rg2;
      Vector3D fc;
      
      #ifdef SFML
      std::cout << "TEST" << std::endl;
      //SFML
      sf::ContextSettings settings;
      settings.depthBits = 24;
      settings.stencilBits = 8;
      settings.antialiasingLevel = 0;
      settings.majorVersion = 0;
      settings.minorVersion = 0;

      sf::Window app(sf::VideoMode(800, 600), "Window", sf::Style::Default, settings);
      sf::Clock clock;
      sf::Event event;
      glInit();

      // input reading depends on how quickly the simulation is running

      // framerate only affects how often the polymer is drawn to the screen the
      // simulation continues even if a frame is not drawn
       TYPE_FLOAT framerate = 24;
      TYPE_FLOAT invFramerate = 1/framerate;
      TYPE_FLOAT lastFrame;
      TYPE_FLOAT scale = 5.5;
      TYPE_FLOAT xCamera = 0;
      TYPE_FLOAT yCamera = 0;
      TYPE_FLOAT xAngle = 0;
      TYPE_FLOAT yAngle = 0;
      TYPE_FLOAT moveSpeed = 0.05;
      TYPE_FLOAT rotateSpeed = 1;

      std::cout<<time(NULL)<<std::endl;
      int iRotateCamera = 0;
      while(w.getStep() < WorldSettings::steps && isFinished != 1)
        {
	  
	  while (app.pollEvent(event))
            {
	      if(event.type == sf::Event::MouseWheelMoved)
                {
		  scale -= event.mouseWheel.delta/4.0;
                }
	      // the window becomes not responding unless the event stack is emptied
            }
	    
	  if(sf::Keyboard::isKeyPressed(sf::Keyboard::Left))
            {
	      xCamera -= moveSpeed;
            }
	  if(sf::Keyboard::isKeyPressed(sf::Keyboard::Right))
            {
	      xCamera += moveSpeed;
            }

	  if(sf::Keyboard::isKeyPressed(sf::Keyboard::Up))
            {
	      yCamera += moveSpeed;
            }
	  if(sf::Keyboard::isKeyPressed(sf::Keyboard::Down))
            {
	      yCamera -= moveSpeed;
            }
	  if(sf::Keyboard::isKeyPressed(sf::Keyboard::A))
            {
	      yAngle -= rotateSpeed;
            }
	  if(sf::Keyboard::isKeyPressed(sf::Keyboard::D))
            {
	      yAngle += rotateSpeed;
            }

	  if(sf::Keyboard::isKeyPressed(sf::Keyboard::W))
            {
	      xAngle -= rotateSpeed;
            }
	  if(sf::Keyboard::isKeyPressed(sf::Keyboard::S))
            {
	      xAngle += rotateSpeed;
	    }
	  if(sf::Keyboard::isKeyPressed(sf::Keyboard::Q))
	    {
	      WorldSettings::terminate = true;
	    }

	  w.simulate();

	  //std::cout<<w.rkInverse()<<std::endl;
	  	  
	  if(WorldSettings::terminate)
	    {
	      w.savePositions();
	      w.saveStateBinary(WorldSettings::checkpointFile);
	      isFinished = 1;
	    }
	  if(WorldSettings::saveCheckpoint)
	    {
	      w.savePositions();
	      w.saveStateBinary(WorldSettings::checkpointFile);
	      WorldSettings::saveCheckpoint = false;
	    }

	  lastFrame = clock.getElapsedTime().asSeconds();
	  if(lastFrame > invFramerate)
            {
	      clock.restart();
	      
	      //draw
	      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	      glMatrixMode(GL_MODELVIEW);
	      glLoadIdentity();

	      // move back so that the entire polymer fits on the screen
	      glTranslatef(xCamera, yCamera, -scale);

	      glRotatef(xAngle, 1, 0, 0);
	      glRotatef(yAngle, 0, 1, 0);

	      //rotate camera around the polymer
	      //glRotatef(iRotateCamera++ /lastFrame* 0.003, 0.f, 1.f, 1.f);

	      //draw axes
	        glColor3f(1.0f,1.0f,1.0f);
                glBegin(GL_LINES);
		glVertex3f(0.f,0.f,0.f);
		glVertex3f(100.f,0.f,0.f);

		glVertex3f(0.f,0.f,0.f);
		glVertex3f(0.f,100.f,0.f);

		glVertex3f(0.f,0.f,0.f);
		glVertex3f(0.f,0.f,100.f);
		glEnd();
	      //
	      w.draw(scale);

	      app.display();
            }
        }
      #endif
    } 

  //clean up
  
  //w.saveContactMap(WorldSettings::contactMapFile);
  /*  if the simulation finishes all steps, it did not receive a SIGTERM,
   *  and the terminate flag was not set, save a checkpoint of the state
   *  at the last file write (will not work with SCREEN)
   */
  if(!WorldSettings::terminate){
    w.saveSeparationIJ(WorldSettings::separationFile);
    w.savePositions();
    w.saveStateBinary(WorldSettings::checkpointFile);
  }
  delete WorldSettings::nd;
  
  for(int i = 0; i < numLines; ++i) {
    delete [] outputArray[i];
  }
  delete [] outputArray;
  
  std::cout<<"LJ Error Count: "<<WorldSettings::errorCountLJ<<" out of "<<WorldSettings::totalCountLJ<<std::endl;
  std::cout<<"FENE Error Count: "<<WorldSettings::errorCountFENE<<" out of "<<WorldSettings::totalCountFENE<<std::endl;
  return isFinished;
}

#ifdef GL
void glInit()
{
  //OpenGL
  glClearDepth(1.f);
  glClearColor(0.f, 0.f, 0.f, 0.f);
  glEnable(GL_DEPTH_TEST);
  //glDepthMask(GL_TRUE);
  glDepthFunc(GL_LEQUAL);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  //Enable Lighting
  //glEnable(GL_LIGHTING);
  //glEnable(GL_LIGHT0);
  //glEnable(GL_NORMALIZE);
  //GLfloat lightpos[] = {1., 1., 1., 0.};
  //glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
  //double amb = .3;
  //GLfloat global_ambient[] = {amb,amb,amb, 0.1};
  //glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
  
  //120 is the viewing angle
  //10000 is the max draw distance
  //#ifdef SFML
  gluPerspective(120.f, 1.33f, 1.f, 10000.f);
  //#endif
}
#endif
