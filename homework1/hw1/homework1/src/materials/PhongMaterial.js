class PhongMaterial extends Material {

    constructor(color, specular, lights, vertexShader, fragmentShader) {
        // let lightMVP = light.CalcLightMVP(translate, scale);
        let lightIntensity1 = lights[0].entity.mat.GetIntensity();
        let lightIntensity2 = [1.0,1.0,1.0];
        if(lights.length>1){
            lightIntensity2 = lights[1].entity.mat.GetIntensity();
        }

        let lightfbo1 = lights[0].entity.fbo;
        let lightfbo2 ;
        if(lights.length>1){
            lightfbo2 = lights[1].entity.fbo;
        }

        super({
            // Phong
            'uSampler': { type: 'texture', value: color },
            'uKs': { type: '3fv', value: specular },
            'uLightIntensity1': { type: '3fv', value: lightIntensity1 },
            // Shadow
            'uShadowMap1': { type: 'texture', value: lightfbo1 },

            'uLightIntensity2': { type: '3fv', value: lightIntensity2 },
            // Shadow
            'uShadowMap2': { type: 'texture', value: lightfbo2 },

            // 'uLightMVP': { type: 'matrix4fv', value: lightMVP },

        }, [], vertexShader, fragmentShader);
    }
}

async function buildPhongMaterial(color, specular, lights, vertexPath, fragmentPath) {


    let vertexShader = await getShaderString(vertexPath);
    let fragmentShader = await getShaderString(fragmentPath);

    return new PhongMaterial(color, specular, lights, vertexShader, fragmentShader);

}