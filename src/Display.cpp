#include "Display.h"

#ifdef SFML

Display::Display(){
  inverseDrawFramerate = 1/24;
  inverseEventFramerate = 1/60;
  cameraZoom = 5.5;
  xCameraPosition = 0;
  yCameraPosition = 0;
  xCameraAngle = 0;
  yCameraAngle = 0;
  cameraMoveSpeed = 0.05;
  cameraRotateSpeed = 1;
  cameraZoomSpeed = 0.25;

  
  sf::ContextSettings settings;
  settings.depthBits = 24;
  settings.stencilBits = 8;
  settings.antialiasingLevel = 0;
  settings.majorVersion = 0;
  settings.minorVersion = 0;
  
  window.create(sf::VideoMode(800, 600), "Window", sf::Style::Default, settings);  
  initGL();

  eventClock.restart();
  drawClock.restart();
}

void Display::handleEvents(){
  if(eventClock.getElapsedTime().asSeconds() < inverseEventFramerate)
    return;
  
  eventClock.restart();
  while (window.pollEvent(event))
    {
      if(event.type == sf::Event::MouseWheelMoved)
	{
	  cameraZoom -= cameraZoomSpeed * event.mouseWheel.delta;
	}
      // the window doesn't respond unless the event stack is emptied
    }
  
  if(sf::Keyboard::isKeyPressed(sf::Keyboard::Left))
    {
      xCameraPosition -= cameraMoveSpeed;
    }
  if(sf::Keyboard::isKeyPressed(sf::Keyboard::Right))
    {
      xCameraPosition += cameraMoveSpeed;
    }
  
  if(sf::Keyboard::isKeyPressed(sf::Keyboard::Up))
    {
      yCameraPosition += cameraMoveSpeed;
    }
  if(sf::Keyboard::isKeyPressed(sf::Keyboard::Down))
    {
      yCameraPosition -= cameraMoveSpeed;
    }
  if(sf::Keyboard::isKeyPressed(sf::Keyboard::A))
    {
      yCameraAngle -= cameraRotateSpeed;
    }
  if(sf::Keyboard::isKeyPressed(sf::Keyboard::D))
    {
      yCameraAngle += cameraRotateSpeed;
    }
  
  if(sf::Keyboard::isKeyPressed(sf::Keyboard::W))
    {
      xCameraAngle -= cameraRotateSpeed;
    }
  if(sf::Keyboard::isKeyPressed(sf::Keyboard::S))
    {
      xCameraAngle += cameraRotateSpeed;
    }
  if(sf::Keyboard::isKeyPressed(sf::Keyboard::Q))
    {
      WorldSettings::terminate = true;
    }
}

void Display::draw(const World& world){
  if(drawClock.getElapsedTime().asSeconds() < inverseDrawFramerate)
    return;
  drawClock.restart();
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glTranslatef(xCameraPosition, yCameraPosition, -cameraZoom);
  glRotatef(xCameraAngle, 1, 0, 0);
  glRotatef(yCameraAngle, 0, 1, 0);
  
  drawAxes();

  glPushMatrix();
  glScalef(1.0/cameraZoom,
	   1.0/cameraZoom,
	   1.0/cameraZoom);
  world.draw();
  glPopMatrix();
  
  window.display();
}

void initGL()
{
  //glClearDepth(1.0);
  //glClearColor(0.f, 0.f, 0.f, 0.f);
  glEnable(GL_DEPTH_TEST);
  //glDepthMask(GL_TRUE);
  //glDepthFunc(GL_LEQUAL);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  //Enable Lighting with defaults
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  //glEnable(GL_NORMALIZE);
  GLfloat lightDiffuse[4] = {0.3f, 0.3f, 0.3f, 1.0f};
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
  
  gluPerspective(120.f, // FOV
		 1.33f, // aspect ratio
		 0.1f, // min draw distance
		 10000.f); // max draw distance
}

void drawAxes(){
  glColor3f(1.0f,1.0f,1.0f);

  glBegin(GL_LINES);
  
  glVertex3f(0.f,0.f,0.f);
  glVertex3f(100.f,0.f,0.f);
  
  glVertex3f(0.f,0.f,0.f);
  glVertex3f(0.f,100.f,0.f);
  
  glVertex3f(0.f,0.f,0.f);
  glVertex3f(0.f,0.f,100.f);
  
  glEnd();
}
#endif
