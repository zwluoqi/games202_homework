attribute vec3 aVertexPosition;
attribute vec3 aNormalPosition;
attribute vec2 aTextureCoord;

uniform mat4 uModelMatrix;
uniform mat4 uViewMatrix;
uniform mat4 uProjectionMatrix;
uniform mat4 uLightMVP; 


varying mat4 vWorldToScreen;
varying mat4 vViewMatrix;
varying mat4 vProjectMatrix;
varying highp vec4 vPosWorld;


void main(void) {

  vPosWorld = uModelMatrix * vec4(aVertexPosition, 1.0);
  vPosWorld = vPosWorld.xyzw / vPosWorld.w;
  vWorldToScreen = uProjectionMatrix * uViewMatrix;
  vViewMatrix = uViewMatrix;
  vProjectMatrix = uProjectionMatrix;

  gl_Position = uProjectionMatrix * uViewMatrix * uModelMatrix * vec4(aVertexPosition, 1.0);
}