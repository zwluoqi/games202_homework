class DirectionalLight {

    constructor(lightIntensity, lightColor, lightPos, focalPoint, lightUp, hasShadowMap, gl) {
        this.mesh = Mesh.cube(setTransform(0, 0, 0, 0.2, 0.2, 0.2, 0));
        this.mat = new EmissiveMaterial(lightIntensity, lightColor);
        this.lightPos = lightPos;
        this.focalPoint = focalPoint;
        this.lightUp = lightUp;

        this.hasShadowMap = hasShadowMap;
        this.fbo = new FBO(gl);
        if (!this.fbo) {
            console.log("无法设置帧缓冲区对象");
            return;
        }
    }

    CalcLightMVP(translate, scale) {
        let lightMVP = mat4.create();
        let modelMatrix = mat4.create();
        let viewMatrix = mat4.create();
        let viewMatrix2 = mat4.create();
        let projectionMatrix = mat4.create();
        let projectionMatrix2 = mat4.create();

        // Model transform
        mat4.identity(modelMatrix);
		mat4.translate(modelMatrix, modelMatrix, translate);
		mat4.scale(modelMatrix, modelMatrix, scale);
        // View transform
        // camera.updateMatrixWorld();
		// mat4.invert(viewMatrix, camera.matrixWorld.elements);
        // gl.gluPerspective(fovl,fovn);
        viewMatrix = this.mat4_lookat(this.fromValues ( this.lightPos), this.fromValues (this.focalPoint),this.fromValues ( this.lightUp));
        mat4.lookAt(viewMatrix2, this.lightPos, this.focalPoint, this.lightUp);
        // Projection transform
        // mat4.copy(projectionMatrix, camera.projectionMatrix.elements);
        let left = -100
        let right = 100
        let bottom = -100
        let top = 100
        let near = -420
        let far = -40
        projectionMatrix = this.mat4_ortho(left,right,bottom,top,near,far)
        mat4.ortho(projectionMatrix2,left,right,bottom,top,near,far);

        mat4.multiply(lightMVP, projectionMatrix, viewMatrix);
        mat4.multiply(lightMVP, lightMVP, modelMatrix);

        return lightMVP;
    }


    fromValues(pos){
        return  vec3.fromValues (pos[0],pos[1],pos[2])

    }
        

    /*
    * eye: the position of the eye point
    * target: the position of the target point
    * up: the direction of the up vector
    *
    * x_axis.x  x_axis.y  x_axis.z  -dot(x_axis,eye)
    * y_axis.x  y_axis.y  y_axis.z  -dot(y_axis,eye)
    * z_axis.x  z_axis.y  z_axis.z  -dot(z_axis,eye)
    *        0         0         0                 1
    *
    * z_axis: normalize(eye-target), the backward vector
    * x_axis: normalize(cross(up,z_axis)), the right vector
    * y_axis: cross(z_axis,x_axis), the up vector
    *
    * see http://www.songho.ca/opengl/gl_camera.html
    */
mat4_lookat(light_eye,light_target,light_up){
    let m = mat4.create();
    mat4.identity(m);

    let subtract = vec3.create() ;
    let tmp = vec3.create() ;
    vec3.subtract(subtract,light_eye , light_target)
    let z = vec3.create() ;
    vec3.normalize(z,subtract);
    let x =vec3.create() ;
    vec3.cross(tmp,light_up, z);
     vec3.normalize(x,tmp);
    let y =vec3.create() ;
    vec3.cross(tmp,z, x);
     vec3.normalize(y,tmp);

    m[0] = x[0];
    m[4] = x[1];
    m[8] = x[2];

    m[1] = y[0];
    m[5] = y[1];
    m[9] = y[2];

    m[2] = z[0];
    m[6] = z[1];
    m[10] = z[2];
    let tmpfloat = 0;
    tmpfloat = vec3.dot(x, light_eye);
    m[12] = -tmpfloat; //相当于原来要位移的，在新的坐标系下是位移多少，有个改变
    tmpfloat = vec3.dot(y, light_eye);
    m[13] = -tmpfloat;
    let tmpfloat2 = z[0]*light_eye[0]+z[1]*light_eye[1]+z[2]*light_eye[2];
    tmpfloat = vec3.dot(z, light_eye);
    m[14] = -tmpfloat;


    return m;
}
// @param {Number} m03 Component in column 0, row 3 position (index 3)

    
/*
* left, right: the coordinates for the left and right clipping planes
* bottom, top: the coordinates for the bottom and top clipping planes
* near, far: the coordinates for the near and far clipping planes
*
* 2/(r-l)        0         0  -(r+l)/(r-l)
*       0  2/(t-b)         0  -(t+b)/(t-b)
*       0        0   2/(n-f)  -(f+n)/(n-f)
*       0        0         0             1
* see http://docs.gl/gl2/glOrtho
*
* note: opgenl map the near plane to -1, far plane to 1,
*       but i map the near plane to 1, far plane to -1.
*       And if near and far is positive it means the plane is behind viewer.
*/
mat4_ortho(left, right, bottom, top,
    near, far)
{
    let x_range = right - left;
    let y_range = top - bottom;
    let z_range = near - far;  //care the different 
    let m = mat4.create();
    mat4.identity(m);

    m[0] = 2 / x_range;
    m[5] = 2 / y_range;
    m[10] = 2 / z_range;
    m[12] = -(left + right) / x_range;
    m[13] = -(bottom + top) / y_range;
    m[14] = -(near + far) / z_range;
    return m;
}

}
