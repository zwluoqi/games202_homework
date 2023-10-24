1.直接光照函数Diffuse
2.实现 Screen Space Ray Tracing:RayMarch
3.实现间接光照:
vec3 sampleDir = SampleHemisphereCos(s,pdf);
sampleDir = normalize(TBN*sampleDir);

p1hit = RayMarch(vPosWorld.xyz,sampleDir,p1);        
if (p1hit > 0.0){

p1hitUV = GetScreenCoordinate(p1);
LLoop = EvalDiffuse(sampleDir,wo,uv)/pdf;
LLoop *= EvalDiffuse(uLightDir,-sampleDir,p1hitUV)*EvalDirectionalLight(p1hitUV) ;
L+=LLoop;

}else{
}
// return;