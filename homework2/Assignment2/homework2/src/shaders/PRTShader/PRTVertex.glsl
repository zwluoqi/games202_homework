attribute vec3 aVertexPosition;
attribute vec3 aNormalPosition;
attribute vec2 aTextureCoord;



attribute mat3 aPrecomputeLT;
// attribute mat3 aPrecomputeLT1;
// attribute mat3 aPrecomputeLT2;



uniform mat3  uLightMat3ColorR;
uniform mat3  uLightMat3ColorG;
uniform mat3  uLightMat3ColorB;

// attribute mat3 aPrecomputeLT2;

uniform mat4 uModelMatrix;
uniform mat4 uViewMatrix;
uniform mat4 uProjectionMatrix;

#define PI 3.141592653589793


varying highp vec3 vNormal;
varying highp vec3 vColor;

void main(void){

    vNormal = (uModelMatrix * vec4(aNormalPosition, 0.0)).xyz;


    mat3 light[3];
    light[0] = uLightMat3ColorR;
    light[1] = uLightMat3ColorG;
    light[2] = uLightMat3ColorB;

    vec3 c = vec3(0,0,0);
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            c.r += aPrecomputeLT[i][j] * uLightMat3ColorR[i][j];
            c.g += aPrecomputeLT[i][j] * uLightMat3ColorG[i][j];
            c.b += aPrecomputeLT[i][j] * uLightMat3ColorB[i][j];
        }
    }
    vColor = c;


  gl_Position = uProjectionMatrix * uViewMatrix * uModelMatrix *
                vec4(aVertexPosition, 1.0);

    
}
