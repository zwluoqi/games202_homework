#ifdef GL_ES
precision highp float;
#endif

uniform vec3 uLightDir;
uniform vec3 uCameraPos;
uniform vec3 uLightRadiance;

uniform sampler2D uGDiffuse;
uniform sampler2D uGDepth;
uniform sampler2D uGNormalWorld;
uniform sampler2D uGShadow;
uniform sampler2D uGPosWorld;

uniform vec3 uReslution;

varying mat4 vWorldToScreen;
varying mat4 vViewMatrix;
varying mat4 vProjectMatrix;
varying highp vec4 vPosWorld;


#define M_PI 3.1415926535897932384626433832795
#define TWO_PI 6.283185307
#define INV_PI 0.31830988618
#define INV_TWO_PI 0.15915494309

#define MAX_ITERATIONS = 20.0,
#define MAX_STRIDE = 16.0,
#define SEARCH_STEPS = 3.0;
#define _MaximumMarchDistance = 1000;

float rand(vec2 co){
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}

float Rand1(inout float p) {
  p = fract(p * .1031);
  p *= p + 33.33;
  p *= p + p;
  return fract(p);
}

vec2 Rand2(inout float p) {
  return vec2(Rand1(p), Rand1(p));
}

vec3 Rand3(inout float p) {
  return vec3(Rand1(p), Rand1(p),Rand1(p));
}

float InitRand(vec2 uv) {
	vec3 p3  = fract(vec3(uv.xyx) * .1031);
  p3 += dot(p3, p3.yzx + 33.33);
  return fract((p3.x + p3.y) * p3.z);
}

vec3 SampleHemisphereUniform(inout float s, out float pdf) {
  vec2 uv = Rand2(s);
  float z = uv.x;
  float phi = uv.y * TWO_PI;
  float sinTheta = sqrt(1.0 - z*z);
  vec3 dir = vec3(sinTheta * cos(phi), sinTheta * sin(phi), z);
  pdf = INV_TWO_PI;
  return dir;
}

vec3 SampleHemisphereCos(inout float s, out float pdf) {
  vec2 uv = Rand2(s);
  float z = sqrt(1.0 - uv.x);
  float phi = uv.y * TWO_PI;
  float sinTheta = sqrt(uv.x);
  vec3 dir = vec3(sinTheta * cos(phi), sinTheta * sin(phi), z);
  pdf = z * INV_PI;
  return dir;
}

void LocalBasis(vec3 n, out vec3 b1, out vec3 b2) {
  float sign_ = sign(n.z);
  if (n.z == 0.0) {
    sign_ = 1.0;
  }
  float a = -1.0 / (sign_ + n.z);
  float b = n.x * n.y * a;
  b1 = vec3(1.0 + sign_ * n.x * n.x * a, sign_ * b, -sign_ * n.x);
  b2 = vec3(b, sign_ + n.y * n.y * a, -n.y);
}

vec4 Project(vec4 a) {
  return a / a.w;
}

float GetDepth(vec3 posWorld) {
  float depth = (vWorldToScreen * vec4(posWorld, 1.0)).w;
  return depth;
}

/*
 * Transform point from world space to screen space([0, 1] x [0, 1])
 *
*/
vec2 GetScreenCoordinate(vec3 posWorld) {
  vec2 uv = Project(vWorldToScreen * vec4(posWorld, 1.0)).xy * 0.5 + 0.5;
  return uv;
}

//[0-1000]
float GetGBufferDepth(vec2 uv) {
  float depth = texture2D(uGDepth, uv).x;
  if (depth < 1e-2) {
    depth = 1000.0;
  }
  return depth;
}

vec3 GetGBufferNormalWorld(vec2 uv) {
  vec3 normal = texture2D(uGNormalWorld, uv).xyz;
  return normal;
}

vec3 GetGBufferPosWorld(vec2 uv) {
  vec3 posWorld = texture2D(uGPosWorld, uv).xyz;
  return posWorld;
}

float GetGBufferuShadow(vec2 uv) {
  float visibility = texture2D(uGShadow, uv).x;
  return visibility;
}

vec3 GetGBufferDiffuse(vec2 uv) {
  vec3 diffuse = texture2D(uGDiffuse, uv).xyz;
  diffuse = pow(diffuse, vec3(2.2));
  return diffuse;
}

mat4 inverse_mat4(mat4 m)
{
    float Coef00 = m[2][2] * m[3][3] - m[3][2] * m[2][3];
    float Coef02 = m[1][2] * m[3][3] - m[3][2] * m[1][3];
    float Coef03 = m[1][2] * m[2][3] - m[2][2] * m[1][3];
    
    float Coef04 = m[2][1] * m[3][3] - m[3][1] * m[2][3];
    float Coef06 = m[1][1] * m[3][3] - m[3][1] * m[1][3];
    float Coef07 = m[1][1] * m[2][3] - m[2][1] * m[1][3];
    
    float Coef08 = m[2][1] * m[3][2] - m[3][1] * m[2][2];
    float Coef10 = m[1][1] * m[3][2] - m[3][1] * m[1][2];
    float Coef11 = m[1][1] * m[2][2] - m[2][1] * m[1][2];
    
    float Coef12 = m[2][0] * m[3][3] - m[3][0] * m[2][3];
    float Coef14 = m[1][0] * m[3][3] - m[3][0] * m[1][3];
    float Coef15 = m[1][0] * m[2][3] - m[2][0] * m[1][3];
    
    float Coef16 = m[2][0] * m[3][2] - m[3][0] * m[2][2];
    float Coef18 = m[1][0] * m[3][2] - m[3][0] * m[1][2];
    float Coef19 = m[1][0] * m[2][2] - m[2][0] * m[1][2];
    
    float Coef20 = m[2][0] * m[3][1] - m[3][0] * m[2][1];
    float Coef22 = m[1][0] * m[3][1] - m[3][0] * m[1][1];
    float Coef23 = m[1][0] * m[2][1] - m[2][0] * m[1][1];
    
    const vec4 SignA = vec4( 1.0, -1.0,  1.0, -1.0);
    const vec4 SignB = vec4(-1.0,  1.0, -1.0,  1.0);
    
    vec4 Fac0 = vec4(Coef00, Coef00, Coef02, Coef03);
    vec4 Fac1 = vec4(Coef04, Coef04, Coef06, Coef07);
    vec4 Fac2 = vec4(Coef08, Coef08, Coef10, Coef11);
    vec4 Fac3 = vec4(Coef12, Coef12, Coef14, Coef15);
    vec4 Fac4 = vec4(Coef16, Coef16, Coef18, Coef19);
    vec4 Fac5 = vec4(Coef20, Coef20, Coef22, Coef23);
    
    vec4 Vec0 = vec4(m[1][0], m[0][0], m[0][0], m[0][0]);
    vec4 Vec1 = vec4(m[1][1], m[0][1], m[0][1], m[0][1]);
    vec4 Vec2 = vec4(m[1][2], m[0][2], m[0][2], m[0][2]);
    vec4 Vec3 = vec4(m[1][3], m[0][3], m[0][3], m[0][3]);
    
    vec4 Inv0 = SignA * (Vec1 * Fac0 - Vec2 * Fac1 + Vec3 * Fac2);
    vec4 Inv1 = SignB * (Vec0 * Fac0 - Vec2 * Fac3 + Vec3 * Fac4);
    vec4 Inv2 = SignA * (Vec0 * Fac1 - Vec1 * Fac3 + Vec3 * Fac5);
    vec4 Inv3 = SignB * (Vec0 * Fac2 - Vec1 * Fac4 + Vec2 * Fac5);
    
    mat4 Inverse = mat4(Inv0, Inv1, Inv2, Inv3);
    
    vec4 Row0 = vec4(Inverse[0][0], Inverse[1][0], Inverse[2][0], Inverse[3][0]);
    
    float Determinant = dot(m[0], Row0);
    
    Inverse /= Determinant;
    
    return Inverse;
}


/*
 * Evaluate diffuse bsdf value.
 *
 * wi, wo are all in world space.
 * uv is in screen space, [0, 1] x [0, 1].
 *
 */
vec3 EvalDiffuse(vec3 wi, vec3 wo, vec2 uv) {
  //经过计算后在世界坐标系下的法线
  vec3 N = GetGBufferNormalWorld(uv);
  
  //定义在mesh文件中的漫反射率
  vec3 diffuse =  GetGBufferDiffuse(uv);

  //vec3 L =  diffuse * max(0.0,dot(N,wo));
  // vec3 L = diffuse * max(0.0,dot(N,wo))*1.0*INV_PI;
  vec3 L = diffuse *INV_PI;
  //vec3 L = vec3(0.0);
  return L;
}

vec3 Diffuse(vec3 wi, vec3 wo, vec2 uv){
 //经过计算后在世界坐标系下的法线
  vec3 N = GetGBufferNormalWorld(uv);
  
  //定义在mesh文件中的漫反射率
  vec3 diffuse =  GetGBufferDiffuse(uv);

  //vec3 L =  diffuse * max(0.0,dot(N,wo));
  vec3 L = diffuse * max(0.0,dot(N,wi))*1.0;
  //vec3 L = vec3(0.0);
  return L;
}

/*
 * Evaluate directional light with shadow map
 * uv is in screen space, [0, 1] x [0, 1].
 *
 */
vec3 EvalDirectionalLight(vec2 uv) {


  float visibility = GetGBufferuShadow(uv) ;
  vec3 Le = visibility*uLightRadiance*1.0;
  //vec3 Le = vec3(0.0);
  return Le;
}



//http://casual-effects.blogspot.com/2014/08/screen-space-ray-tracing.html

float RayMarch(vec3 ori, vec3 dir, out vec3 hitPos) {
    vec2 screenSize = uReslution.xy;//2000*1000
    float uWidth = screenSize.x;
    float uHeight = screenSize.y;
    const float nearPlaneZ = -0.001;
    const float maxSteps = 1000.0;
    const float maxRayTraceDistance = 1000.0;
    float stride = 1.0 * max(1.0/uWidth, 1.0/uHeight);
    const float jitterFraction = 1.0;
    const float csZThickness = 0.1;

    // Transform into camera space 
    vec4 tmpVec4 = vViewMatrix * vec4(ori, 1.0);
    vec3 csOrigin = tmpVec4.xyz/tmpVec4.w;
    tmpVec4 = vViewMatrix * vec4(dir, 0.0);
    vec3 csDirection = normalize(tmpVec4.xyz);

    // Clip ray to a near plane in 3D
    float rayLength = ((csOrigin.z + csDirection.z * maxRayTraceDistance) > nearPlaneZ) ? 
                    (nearPlaneZ - csOrigin.z) / csDirection.z : maxRayTraceDistance;
    vec3 csEndPoint = csOrigin + csDirection * rayLength;

    // Project into screen space
    vec4 H0 = vProjectMatrix * vec4(csOrigin, 1.0);
    vec4 H1 = vProjectMatrix * vec4(csEndPoint, 1.0);
    float k0 = 1.0 / H0.w;
    float k1 = 1.0 / H1.w;

    // Switch the original points to values that interpolate linearly in 2D
    vec3 Q0 = csOrigin * k0; 
    vec3 Q1 = csEndPoint * k1;

    // Screen-space endpoints
    vec2 P0 = vec2(H0.xy * k0) * 0.5 + 0.5;
    vec2 P1 = vec2(H1.xy * k1) * 0.5 + 0.5;
    // vec2 P0 = vec2(H0.xy * k0) ;
    // vec2 P1 = vec2(H1.xy * k1) ;

    vec2 hitPixel = vec2(-1.0, -1.0);

    // If the line is degenerate, make it cover at least one pixel
    // to avoid handling zero-pixel extent as a special case later
    P1 += vec2((dot(P0 - P1, P0 - P1) < 0.0001) ? 0.01 : 0.0);

    vec2 delta = P1 - P0;

    // Permute so that the primary iteration is in x to reduce
    // large branches later
    bool permute = false;
	if (abs(delta.x) < abs(delta.y)) 
    {
		// More-vertical line. Create a permutation that swaps x and y in the output
		permute = true;
        // Directly swizzle the inputs
		delta = delta.yx;
		P1 = P1.yx;
		P0 = P0.yx;        
	}

    // From now on, "x" is the primary iteration direction and "y" is the secondary one
    float stepDirection = sign(delta.x);
    float invdx = stepDirection / delta.x;
    vec2 dP = vec2(stepDirection, invdx * delta.y);

    // Track the derivatives of Q and k
    vec3 dQ = (Q1 - Q0) * invdx;
    float dk = (k1 - k0) * invdx;

    // Scale derivatives by the desired pixel stride
	dP *= stride; 
    dQ *= stride; 
    dk *= stride;

    // Offset the starting values by the jitter fraction
	P0 += dP * jitterFraction; 
    Q0 += dQ * jitterFraction; 
    k0 += dk * jitterFraction;

    // Slide P from P0 to P1, (now-homogeneous) Q from Q0 to Q1, and k from k0 to k1
    vec3 Q = Q0;
    float k = k0;

	// We track the ray depth at +/- 1/2 pixel to treat pixels as clip-space solid 
	// voxels. Because the depth at -1/2 for a given pixel will be the same as at 
	// +1/2 for the previous iteration, we actually only have to compute one value 
	// per iteration.
	float prevZMaxEstimate = csOrigin.z;
    float stepCount = 0.0;
    float rayZMax = prevZMaxEstimate, rayZMin = prevZMaxEstimate;
    float sceneZMax = csOrigin.z + 1e4;

    // P1.x is never modified after this point, so pre-scale it by 
    // the step direction for a signed comparison
    float end = P1.x * stepDirection;

    // We only advance the z field of Q in the inner loop, since
    // Q.xy is never used until after the loop terminates.
    vec2 P = P0;
    for(int i = 0;i < int(maxSteps); ++i)
    {
        if(P.x < 0.0 || P.x > 1.0 || P.y < 0.0 || P.y > 1.0)
        {
            break;
        }

        if(!(((P.x * stepDirection) <= end) && 
             ((rayZMax < sceneZMax - csZThickness) ||
              (rayZMin > sceneZMax)) &&
             (sceneZMax != 0.0)))
        {
            break;
        }

        hitPixel = permute ? P.yx : P;
        // The depth range that the ray covers within this loop
        // iteration.  Assume that the ray is moving in increasing z
        // and swap if backwards.  Because one end of the interval is
        // shared between adjacent iterations, we track the previous
        // value and then swap as needed to ensure correct ordering
        rayZMin = prevZMaxEstimate;

        // Compute the value at 1/2 pixel into the future
        rayZMax = (dQ.z * 0.5 + Q.z) / (dk * 0.5 + k);
        prevZMaxEstimate = rayZMax;

        if (rayZMin > rayZMax) 
        { 
            // Swap
            float tmp = rayZMin;
            rayZMin = rayZMax;
            rayZMax = tmp;
        }

        // Camera-space z of the background
        sceneZMax = -texture2D(uGDepth, hitPixel).r;

        P += dP, Q.z += dQ.z, k += dk, stepCount += 1.0;
    }

    Q.xy += dQ.xy * stepCount;
    hitPos = Q * (1.0 / k);

    tmpVec4 = (inverse_mat4(vViewMatrix)*vec4(hitPos,1.0));
    hitPos = tmpVec4.xyz/tmpVec4.w;

    return (rayZMax >= sceneZMax - csZThickness) && (rayZMin <= sceneZMax)?1.0:0.0;
  
}


float RayMarch2(vec3 ori, vec3 dir, out vec3 hitPos) {
  hitPos = vec3(0.0);

  float maxDistance = 1000.0;
  float resolution  = 0.3;
  int   steps       = 200;
  float thickness   = 1.0;
  
  vec2 texSize  = uReslution.xy;//2000*1000

  vec3 end = ori + dir*maxDistance;

  vec3 nextPos = ori;

  vec3  increment = dir*0.1;
  
  float search0 = 0.0;
  float search1 = 0.0;

  int hit0 = 0;
  int hit1 = 0;

  float depth        = thickness;
  float gbuffDepth = thickness;
  float posDepth = thickness;

   for (int i  = 0; i < 100; ++i) {

     nextPos += increment;
      vec2 uv = GetScreenCoordinate(nextPos.xyz);

    gbuffDepth = GetGBufferDepth(uv);
    posDepth = GetDepth(nextPos.xyz);
    depth = posDepth - gbuffDepth;

    if (depth > 1e-2 && depth < thickness) {
      hit0 = 1;
      break;
    } else {
      search0 = search1;
    }
   }
  hitPos = nextPos;
  return float(hit0);
  
}

#define SAMPLE_NUM 1

void main() {

  vec3 L = vec3(0.0);
  // L = GetGBufferDiffuse(GetScreenCoordinate(vPosWorld.xyz));
  // vec3 color = pow(clamp(L, vec3(0.0), vec3(1.0)), vec3(1.0 / 2.2));
  // gl_FragColor = vec4(vec3(color.rgb), 1.0);

  vec2 uv = GetScreenCoordinate(vPosWorld.xyz);
  float s = InitRand(gl_FragCoord.xy);

  vec3 wo = normalize(uCameraPos.xyz - vPosWorld.xyz );
  vec3 color = Diffuse(uLightDir,wo,uv);
  vec3 shadow = EvalDirectionalLight(uv);
  color *= shadow;
  // float depth = GetGBufferDepth(uv);

  vec3 normal = normalize(GetGBufferNormalWorld(uv));
  // vec3 dir = normalize(reflect(wo, normal));
  // float dN2W = dot(wo,normal);
  // vec3 p0;
  // float p0hit =0.0;
  // vec3 p0hitColor;
  // vec2 p0hitUV;
  // p0hit = RayMarch(vPosWorld.xyz,-dir,p0);
  // if(p0hit > 0.0){
  //   p0hitUV = GetScreenCoordinate(p0);
  //   p0hitColor = Diffuse(dir,wo,p0hitUV)*2.0;
  //   // p0hitColor = abs((p0-p00).xyz)/1.0;
  // }else{

  // }

  vec3 b1,b2;
  LocalBasis(normal,b1,b2);
  mat3 TBN = mat3(b1,b2,normal);

  float pdf = 1.0;
  vec3 LLoop = vec3(0.0);
  float visibility = GetGBufferuShadow(uv) ;
  if(visibility <= 0.3)
  {
      float ignoreCount = 0.0;
      vec3 p1;
      float p1hit =0.0;
      vec2 p1hitUV;
      for(int i=0;i<SAMPLE_NUM;i++){
          vec3 sampleDir = SampleHemisphereCos(s,pdf);
          sampleDir = normalize(TBN*sampleDir);
          
          p1hit = RayMarch(vPosWorld.xyz,sampleDir,p1);        
          if (p1hit > 0.0){
            // float dis = distance(p1,vPosWorld.xyz);
            // float factor = max(0.0, 1.0/sqrt(0.7+dis)-0.2);

            p1hitUV = GetScreenCoordinate(p1);
            LLoop = EvalDiffuse(sampleDir,wo,uv)/pdf;
            LLoop *= Diffuse(uLightDir,-sampleDir,p1hitUV)*EvalDirectionalLight(p1hitUV)*1.0*1.0 ;
            L+=LLoop;

            // vec3 disDir = p1 - vPosWorld.xyz;


            //1/sqrt(x+0.2)-1
            // gl_FragColor = vec4((Diffuse(uLightDir,-sampleDir,p1hitUV)*factor), 1.0);

          }else{
            // ignoreCount++;   
            // gl_FragColor = vec4(0.0,0.0,0.0,1.0);
          }
          // return;
      }
      float counter = float(SAMPLE_NUM-int(ignoreCount));
      if(counter>0.0){
        L = max(L/float(counter),vec3(0.0,0.0,0.0));
      }
   }
 
  //  float visibility = GetGBufferuShadow(uv) ;

  gl_FragColor = vec4(vec3(color.rgb), 1.0);
  // gl_FragColor = vec4(L+vec3(color.rgb), 1.0);
  // gl_FragColor = vec4(p0hitColor, 1.0);
  // gl_FragColor = vec4(vec3(p0hit,0.0,0.0), 1.0);
  // gl_FragColor = vec4(L, 1.0);
  gl_FragColor = vec4(L+vec3(color.rgb), 1.0);

  // gl_FragColor = vec4(vec3(color.rgb)+hitColor*visibale, 1.0);
}
