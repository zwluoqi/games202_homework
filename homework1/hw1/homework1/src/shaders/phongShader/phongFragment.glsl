#ifdef GL_ES
precision mediump float;
#endif

// Phong related variables
uniform sampler2D uSampler;
uniform vec3 uKd;
uniform vec3 uKs;
uniform vec3 uLightPos1;
uniform vec3 uLightPos2;
uniform vec3 uCameraPos;
uniform vec3 uLightIntensity1;
uniform vec3 uLightIntensity2;

varying highp vec2 vTextureCoord;
varying highp vec3 vFragPos;
varying highp vec3 vNormal;

// Shadow map related variables
#define NUM_SAMPLES 20
#define BLOCKER_SEARCH_NUM_SAMPLES NUM_SAMPLES
#define PCF_NUM_SAMPLES NUM_SAMPLES
#define NUM_RINGS 10

#define EPS 1e-3
#define PI 3.141592653589793
#define PI2 6.283185307179586

uniform sampler2D uShadowMap1;
uniform sampler2D uShadowMap2;

varying vec4 vPositionFromLight1;
varying vec4 vPositionFromLight2;

#define LIGHT_WIDTH 2

highp float rand_1to1(highp float x ) { 
  // -1 -1
  return fract(sin(x)*10000.0);
}

highp float rand_2to1(vec2 uv ) { 
  // 0 - 1
	const highp float a = 12.9898, b = 78.233, c = 43758.5453;
	highp float dt = dot( uv.xy, vec2( a,b ) ), sn = mod( dt, PI );
	return fract(sin(sn) * c);
}

float unpack(vec4 rgbaDepth) {
    const vec4 bitShift = vec4(1.0, 1.0/256.0, 1.0/(256.0*256.0), 1.0/(256.0*256.0*256.0));
    return dot(rgbaDepth, bitShift);
}

// float unpack(vec4 rgbaDepth) {
//     // const vec4 bitShift = vec4(1.0, 1.0/256.0, 1.0/(256.0*256.0), 1.0/(256.0*256.0*256.0));
//     // return dot(rgbaDepth, bitShift);
//     return rgbaDepth.r;
// }

vec2 poissonDisk[NUM_SAMPLES];

void poissonDiskSamples( const in vec2 randomSeed ) {

  float ANGLE_STEP = PI2 * float( NUM_RINGS ) / float( NUM_SAMPLES );
  float INV_NUM_SAMPLES = 1.0 / float( NUM_SAMPLES );

  float angle = rand_2to1( randomSeed ) * PI2;
  float radius = INV_NUM_SAMPLES;
  float radiusStep = radius;

  for( int i = 0; i < NUM_SAMPLES; i ++ ) {
    poissonDisk[i] = vec2( cos( angle ), sin( angle ) ) * pow( radius, 0.75 );
    radius += radiusStep;
    angle += ANGLE_STEP;
  }
}

void uniformDiskSamples( const in vec2 randomSeed ) {

  float randNum = rand_2to1(randomSeed);
  float sampleX = rand_1to1( randNum ) ;
  float sampleY = rand_1to1( sampleX ) ;

  float angle = sampleX * PI2;
  float radius = sqrt(sampleY);

  for( int i = 0; i < NUM_SAMPLES; i ++ ) {
    poissonDisk[i] = vec2( radius * cos(angle) , radius * sin(angle)  );

    sampleX = rand_1to1( sampleY ) ;
    sampleY = rand_1to1( sampleX ) ;

    angle = sampleX * PI2;
    radius = sqrt(sampleY);
  }
}

vec2 findBlocker( sampler2D ShadowMap,  vec2 uv, float zReceiver ) {
  // vec3 light2fragPos =  uLightPos - vPositionFromLight.xyz;	
  float avgBlockerDepth = 0.0;
	float numBlockers = 0.0;
	float blockerSum = 0.0;

  //找到周围点中的遮挡该目标点的点，记录其与光源的距离。
  for (int i = 0; i < BLOCKER_SEARCH_NUM_SAMPLES; i++)
	{
      vec2 shadowCoord = poissonDisk[i]*float(BLOCKER_SEARCH_NUM_SAMPLES)/2048.0+uv;
      vec4 rgbaDepth = texture2D(ShadowMap,shadowCoord);
      float depthZ = unpack(rgbaDepth);

      //被挡住了
      if(zReceiver > depthZ && depthZ>0.0){
        blockerSum += depthZ;
        numBlockers += 1.0;
      }
  }

  avgBlockerDepth = blockerSum / numBlockers;

	return vec2(avgBlockerDepth, numBlockers);
}


float useShadowMap(sampler2D ShadowMap, vec4 shadowCoord){
  //return 1.0;
  vec4 rgbaDepth = texture2D(ShadowMap,shadowCoord.xy);
  float depthZ = unpack(rgbaDepth);
  //depthZ = depthZ*0.5+0.5;
  float c = abs(shadowCoord.z -depthZ);
  if(c < 0.01){//biase
    return 1.0;
  }
  float shadow = shadowCoord.z > depthZ  ? 0.0 : 1.0;
  return shadow;
}


float PCF(sampler2D ShadowMap, vec4 coords,float filterRadiusUV) {  
  float invisiableSum =0.0;
  for( int i = 0; i < PCF_NUM_SAMPLES; i ++ ) {
    vec4 shadowCoord = vec4(poissonDisk[i]*filterRadiusUV*float(PCF_NUM_SAMPLES)/2048.0+coords.xy,coords.zw);
    vec4 rgbaDepth = texture2D(ShadowMap,shadowCoord.xy);
    float depthZ = unpack(rgbaDepth);
    //边界,说明光线没有射到物体，可见性为1,避免边界产生阴影,
    if(depthZ<=0.0){
      invisiableSum += 1.0;
    }else{
      invisiableSum +=  useShadowMap(ShadowMap,shadowCoord);
    }
  }
  invisiableSum = invisiableSum/float(PCF_NUM_SAMPLES);
  // if(invisiableSum<0.01){
  //   return 0.0 ;
  // }
  return invisiableSum;
}

float PCSS(sampler2D ShadowMap, vec4 coords){

  // STEP 1: avgblocker depth
  vec2 blockerInfo = findBlocker(ShadowMap,coords.xy,coords.z);
  if (blockerInfo.y < 1.0)
  {
    return 1.0;
  }
  // STEP 2: penumbra size
  float penumbra = abs( coords.z - blockerInfo.x);
  float Wpenumbra = float(LIGHT_WIDTH)*(penumbra)/blockerInfo.x;
  // return Wpenumbra;
  // STEP 3: filtering
  float shadow = PCF(ShadowMap,coords,Wpenumbra);

  return shadow;

}


vec3 blinnPhong(vec3 LightPos,vec3 LightIntensity) {
  vec3 color = texture2D(uSampler, vTextureCoord).rgb;
  color = pow(color, vec3(2.2));

  vec3 ambient = 0.05 * color;

  vec3 lightDir = normalize(LightPos);
  vec3 normal = normalize(vNormal);
  float diff = max(dot(lightDir, normal), 0.0);
  vec3 light_atten_coff =
      LightIntensity / pow(length(LightPos - vFragPos), 2.0);
  vec3 diffuse = diff * light_atten_coff * color;

  vec3 viewDir = normalize(uCameraPos - vFragPos);
  vec3 halfDir = normalize((lightDir + viewDir));
  float spec = pow(max(dot(halfDir, normal), 0.0), 32.0);
  vec3 specular = uKs * light_atten_coff * spec;

  vec3 radiance = (ambient + diffuse + specular);
  vec3 phongColor = pow(radiance, vec3(1.0 / 2.2));
  return phongColor;
}

vec4 GetphongColor(sampler2D ShadowMap,vec4 PositionFromLight,vec3 LightPos,vec3 LightIntensity){

  vec4 shadowCoord = PositionFromLight/PositionFromLight.w;
  shadowCoord.xyz = shadowCoord.xyz*0.5+vec3(0.5,0.5,0.5);
  shadowCoord.xyz = clamp(shadowCoord.xyz,vec3(0.0,0.0,0.0),vec3(1.0,1.0,1.0));
  //vec4 rgbaDepth = texture2D(uShadowMap,shadowCoord.xy);
  //float depthZ = unpack(rgbaDepth);
  //depthZ = depthZ*0.5+0.5;
  //gl_FragColor = vec4(depthZ,0,0,1);
  //gl_FragColor = vec4(shadowCoord.z,0,0,1);
  //return;
  
  float visibility;
  //visibility = useShadowMap(uShadowMap, shadowCoord);
  //poissonDiskSamples(shadowCoord.xy);
  //gl_FragColor = vec4(poissonDisk[18],0,1);
  //return;
  poissonDiskSamples(shadowCoord.xy);
  //visibility = PCF(uShadowMap, shadowCoord,1.0);

  visibility = PCSS(ShadowMap, shadowCoord);
  // gl_FragColor = vec4(visibility,0,0,1);
  // return;
  vec3 phongColor = blinnPhong(LightPos,LightIntensity);
  return vec4(phongColor,visibility);
}

void main(void) {

  // vec4 shadowCoord = vPositionFromLight1/vPositionFromLight1.w;
  // shadowCoord.xyz = shadowCoord.xyz*0.5+vec3(0.5,0.5,0.5);
  // shadowCoord.xyz = clamp(shadowCoord.xyz,vec3(0.0,0.0,0.0),vec3(1.0,1.0,1.0));
  // vec4 rgbaDepth = texture2D(uShadowMap1,shadowCoord.xy);
  // // float depthZ = unpack(rgbaDepth);
  // // depthZ = depthZ*0.5+0.5;
  // gl_FragColor = vec4(rgbaDepth.x*2.0,0,0,1);
  // // gl_FragColor = vec4(shadowCoord.xy,0,1);
  // return;

  vec4 phongColor = GetphongColor(uShadowMap1,vPositionFromLight1,uLightPos1,uLightIntensity1);

  vec4 phongColor2 = GetphongColor(uShadowMap2,vPositionFromLight2,uLightPos2,uLightIntensity2);


  vec4 resultColor = vec4(phongColor.xyz * phongColor.w, 1.0);
  vec4 resultColor2 = vec4(phongColor2.xyz * phongColor2.w, 1.0);

  gl_FragColor = (resultColor + resultColor2);
  // gl_FragColor = resultColor;
  // gl_FragColor = vec4(1.0,1.0,1.0, 1.0);
}