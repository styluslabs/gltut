#include <iostream>

#define GLEW_STATIC
#include <GL/glew.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>

// Eigen is a popular alternative to GLM for more advanced linear algebra
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "glpp.h"

#define QUAD_VERT_SIZE 4
GLfloat quadVertices[] = {
  -1.0f,  1.0f,  0.0f, 1.0f,
   1.0f,  1.0f,  1.0f, 1.0f,
   1.0f, -1.0f,  1.0f, 0.0f,

   1.0f, -1.0f,  1.0f, 0.0f,
  -1.0f, -1.0f,  0.0f, 0.0f,
  -1.0f,  1.0f,  0.0f, 1.0f
};

// changes in newer versions of GLSL:
// - in/out is used instead of varying/attribute
// - the special vars gl_* are not used

const char* toyVert =
R"(#version 120

attribute vec2 position;
//attribute vec2 texcoord;

//varying vec2 Texcoord;

void main()
{
   //Texcoord = texcoord;
   gl_Position = vec4(position, 0.0, 1.0);
}
)";

const char* toyFrag =
R"(#version 120
// Raymarching with fancy lighting
// Refs:
// - http://www.iquilezles.org/www/articles/raymarchingdf/raymarchingdf.htm
// - http://9bitscience.blogspot.com/2013/07/raymarching-distance-fields_14.html
// - https://github.com/nopjia/WebGL-RayMarch
// - http://www.letsdive.in/2014/05/18/glsl---raymarching/
// - http://blog.ruslans.com/2015/01/raymarching-christmas-tree.html

// antialiasing - model ray as a cone (radius r ~ dist from camera).  If object lies within cone (world(p) < r),
//  calculate color and factor into pixel color based on amount of overlap (~(r - world(p))^2)
// - to limit to edges, we should keep track of closest point found within cone as we trace, then factor that in when ray terminates

// geometry operations (union, intersect, etc.)

#define in

uniform vec2 iResolution;
uniform float iGlobalTime;

vec2 opU(in vec2 obj0,in vec2 obj1)
{
  return (obj0.x < obj1.x) ? obj0 : obj1;
}

vec3 opRep(in vec3 p, in vec3 c)
{
  return mod(p,c) - 0.5*c;
}

// scene objects - all are centered at (0,0,0) - world() fn should translate p before calling to translate object position

//Floor
float obj0(in vec3 p){
  return p.y;
}

// IQs RoundBox - see more at http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
float roundBox(in vec3 p, in vec3 dims, in float r) {
  return length(max(abs(p) - dims, 0.0)) - r;
}

float sphere(in vec3 p, in float r) {
  return length(p) - r;
}

// return color based on object index i and position p
vec3 colors(in float i, in vec3 p) {
  if(i == 0.0) {
    // checkerboard floor
    float f = mod(floor(1.0*p.z) + floor(1.0*p.x), 2.0);
    return 0.4 + 0.1*f*vec3(1.0);
  }
  else if(i == 1.0)
    return vec3(1.0,0.5,0.2);
  else if(i == 2.0)
    return vec3(1,0,0);
  else
    return vec3(0.5,0.5,1.0);
}

// scene definition - union of all objects
vec2 world(in vec3 p)
{
  vec2 res = vec2(obj0(p), 0);
  res = opU(res, vec2(roundBox(p - vec3(0,3,0), vec3(1,1,1), 0.25), 1));
  res = opU(res, vec2(roundBox(p - vec3(-2,1.25,4), vec3(1,1,1), 0.25), 2));  //opRep(p, vec3(4,4,4))
  res = opU(res, vec2(sphere(p - vec3(3,3,0), 1.5), 3));
  return res;
}

// end scene definition

// dist < precis is considered a collision
const float precis = 0.01;

// calculate softened shadow or reflection
// march from point p in direction dir (direction to light source for shadows, reflected view direction for relections)
// from http://www.iquilezles.org/www/articles/rmshadows/rmshadows.htm
float softshadow(in vec3 p, in vec3 dir)
{
  const float tmin = 0.1;
  const float tmax = 25.0;
  const float hardness = 16.0;  // larger value -> sharper shadow

  float res = 1.0;
  float t = tmin;
  for(int i = 0; i < 16; i++) {
    float h = world(p + dir*t).x;
    // smaller res means darker shadow
    //  the closer our ray passes to the occluding object (smaller h), the darker the shadow, so res ~ h
    //  dividing by t blurs the edges of the shadows, by allowing a larger h for more distant edges
    res = min(res, hardness*h/t);
    t += clamp(h, 0.2, 1.0);
    if(h < precis || t > tmax)
      break;
  }
  return clamp(res, 0.0, 1.0);
}

// ambient occlusion - 0 (full occlusion) to 1 (none)
// march from point p along normal direction n
float ambientOcclusion(in vec3 p, in vec3 n)
{
  const float range = 0.25;

  float occ = 0.0;
  float sca = 1.0;
  for(int i = 0; i < 5; i++)
  {
    float t = 0.01 + range*float(i)/4.0;
    float h = world(p + n*t).x;
    // if there is a point on another surface closer to the current position than p, so t - h > 0, we add some occlusion
    occ += (t - h)*sca;
    sca *= 0.95;  // weaken effect with distance
  }
  return clamp(1.0 - 3.0*occ, 0.0, 1.0);
}

void main(void)
{
  float time = iGlobalTime;

  // Camera setup
  vec3 cameraPos = vec3(-sin(time)*8.0,4,cos(time)*8.0);  // camera position
  vec3 cameraTarget = vec3(0,0,0);  // point towards which the camera is pointing
  vec3 cameraUp = vec3(0,1,0);  // vector defining up direction of camera view
  const float zoom = 1.0;  // zoom factor - distance to viewing screen used to set ray directions

  // View direction calculation
  // normalized pixel position (-1 to 1)
  vec2 screenXY = -1.0 + 2.0*gl_FragCoord.xy/iResolution.xy;
  // ... adjusted to match aspect ratio
  screenXY *= vec2(iResolution.x/iResolution.y,1.0);

  vec3 cameraDir = normalize(cameraTarget - cameraPos);
  // basis vectors of viewing screen normal to camera direction
  vec3 u = normalize(cross(cameraDir, cameraUp));
  vec3 v = cross(u, cameraDir);
  //vec3 viewDir = normalize(mat3(u, v, cameraDir) * vec3(screenXY, zoom));  // equivalent matrix formulation
  vec3 viewDir = normalize(zoom*cameraDir + screenXY.x*u + screenXY.y*v);  // view ray direction

  // Raymarching
  const vec3 eps = vec3(0.1,0,0);
  const float maxd = 60.0; //Max depth

  vec2 s = vec2(0.1,0.0);
  vec3 p;  // ray end point (collision point after raymarching loop)

  float f = 1.0;
  for(int i = 0; i < 256; i++) {
    if(abs(s.x) < precis || f > maxd)
      break;
    f += s.x;
    p = cameraPos + viewDir*f;
    s = world(p);
  }

  if(f < maxd) {
    // lookup color using index in s.y
    vec3 baseColor = colors(s.y, p);
    // gradient of distance field at p gives surface normal of object
    vec3 n = normalize(vec3(s.x - world(p - eps.xyy).x,
                            s.x - world(p - eps.yxy).x,
                            s.x - world(p - eps.yyx).x));

    // direction *to* light source
    vec3 lightDir = normalize(vec3(0.6, 0.7, 0.5));  //normalize(lightPos - p);
    // reflect view vector
    vec3 viewReflected = reflect(viewDir, n);

    // diffuse light intensity ~ dot(normal, light_direction); light_direction ~ lightPos - p
    float diffuse = clamp(dot(n, lightDir), 0.0, 1.0);

    // shadow
    diffuse *= softshadow(p, lightDir);

    // specular reflection of light source
    float specular = diffuse*pow(clamp(dot(viewReflected, lightDir), 0.0, 1.0), 8.0);

    // ambient occlusion
    float occlusion = ambientOcclusion(p, n);

    // reflections of other objects
    // consider multiplying by smoothstep(-0.1, 0.1, viewReflected.y)
    float reflection = softshadow(p, viewReflected);

    // ambient illumination from above
    float ambient = clamp(0.5 + 0.5*n.y, 0.0, 1.0); //0.25;

    // diffuse lighting from a dim light at the camera position
    float diffuseCam = clamp(dot(n, normalize(cameraPos - p)), 0.0, 1.0);

    // approximate Fresnel effect - brightens edges of objects, where normal is nearly parallel to view direction
    // ref: http://http.developer.nvidia.com/CgTutorial/cg_tutorial_chapter07.html
    float fresnel = pow(clamp(1.0 + dot(n, viewDir), 0.0, 1.0), 6.0);

    // https://www.shadertoy.com/view/Xds3zN includes a term (bac) for light reflected by y = 0 plane

    // fade to black far from camera (disabled)
    float distFade = 1.0; //(1.0 - 0.01*f);
    // light source attenuation
    //float attenuation = 1.0 / (1.0 + pow(length(lightPos - p), 2));

    // combine lighting terms
    float lighting = diffuse; // + specular + occlusion*(0.3*ambient + 0.25*fresnel + 0.1*reflection) + 0.25*diffuseCam;

    vec3 color = baseColor*lighting*distFade;
    //color = pow(color, vec3(1.0/2.2));  // gamma correction
    gl_FragColor = vec4(color, 1.0);
  }
  else
    gl_FragColor = vec4(0.8, 0.9, 1.0, 1.0); // background color
}
)";

int main(int argc, char* argv[])
{
  SDL_Init(SDL_INIT_VIDEO);

  SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);  //3
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);  //2

  SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);

  SDL_Window* window = SDL_CreateWindow("OpenGL", 100, 100, 800, 600, SDL_WINDOW_OPENGL);
  SDL_GLContext context = SDL_GL_CreateContext(window);

  glewExperimental = GL_TRUE;
  glewInit();

  VertexBuffer* vbScreen = new VertexBuffer(quadVertices, sizeof(quadVertices), QUAD_VERT_SIZE);
  VertexArray* vaScreen = new VertexArray;

  ShaderProgram* spScreen = new ShaderProgram(toyVert, toyFrag);
  vaScreen->bindAttrib(spScreen, vbScreen, "position", 2, 0);
  //vaScreen->bindAttrib(spScreen, vbScreen, "texcoord", 2, 2);

  glUniform2f(spScreen->uniform("iResolution"), 800, 600);
  GLuint uniTime = spScreen->uniform("iGlobalTime");

  vaScreen->bind();
  spScreen->use();

  checkGLError();

  SDL_Event windowEvent;
  while(1) {
    if(SDL_PollEvent(&windowEvent)) {
      if (windowEvent.type == SDL_QUIT) break;
      else if(windowEvent.type == SDL_KEYUP && windowEvent.key.keysym.sym == SDLK_ESCAPE) break;
    }

    // draw here
    glUniform1f(uniTime, SDL_GetTicks()/1000.0);  // ms to sec
    glDrawArrays(GL_TRIANGLES, 0, 6);
    // end drawing

    SDL_GL_SwapWindow(window);
  }

  delete spScreen;
  delete vbScreen;
  delete vaScreen;

  // SDL cleanup
  SDL_GL_DeleteContext(context);
  SDL_Quit();
  return 0;
}
