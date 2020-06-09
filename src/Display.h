#ifndef DISPLAY_H
#define DISPLAY_H

#ifdef SFML
#include <SFML/System.hpp>
#include <SFML/Window.hpp>
#include <SFML/OpenGL.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/Config.hpp>

#include <GL/gl.h>
#include <GL/glu.h>

#include "World.h"
#include "WorldSettings.h"

void initGL();
void drawAxes();

class Display{
 public:
  Display();
  void handleEvents();
  void draw(const World& world);
 private:
  sf::Clock eventClock;
  sf::Clock drawClock;
  sf::Event event;
  sf::Window window;

  // framerate only affects how often the polymer is drawn to the screen the
  // simulation continues even if a frame is not drawn
  TYPE_FLOAT inverseDrawFramerate;
  TYPE_FLOAT inverseEventFramerate;
  TYPE_FLOAT cameraZoom;
  TYPE_FLOAT xCameraPosition;
  TYPE_FLOAT yCameraPosition;
  TYPE_FLOAT xCameraAngle;
  TYPE_FLOAT yCameraAngle;
  TYPE_FLOAT cameraMoveSpeed;
  TYPE_FLOAT cameraRotateSpeed;
  TYPE_FLOAT cameraZoomSpeed;

};

#endif

#endif
