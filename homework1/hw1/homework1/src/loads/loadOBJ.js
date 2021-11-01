function loadOBJ(renderer, path, name, objMaterial, transform,autoMove) {

	const manager = new THREE.LoadingManager();
	manager.onProgress = function (item, loaded, total) {
		console.log(item, loaded, total);
	};

	function onProgress(xhr) {
		if (xhr.lengthComputable) {
			const percentComplete = xhr.loaded / xhr.total * 100;
			console.log('model ' + Math.round(percentComplete, 2) + '% downloaded');
		}
	}
	function onError(e ) {
		console.error(e.message)
	 }

	new THREE.MTLLoader(manager)
		.setPath(path)
		.load(name + '.mtl', function (materials) {
			materials.preload();
			new THREE.OBJLoader(manager)
				.setMaterials(materials)
				.setPath(path)
				.load(name + '.obj', function (object) {
					object.traverse(function (child) {
						if (child.isMesh) {
							let geo = child.geometry;
							let mat;
							if (Array.isArray(child.material)) mat = child.material[0];
							else mat = child.material;

							var indices = Array.from({ length: geo.attributes.position.count }, (v, k) => k);
							let mesh = new Mesh({ name: 'aVertexPosition', array: geo.attributes.position.array },
								{ name: 'aNormalPosition', array: geo.attributes.normal.array },
								{ name: 'aTextureCoord', array: geo.attributes.uv.array },
								indices, transform);

							let colorMap = new Texture();
							if (mat.map != null) {
								colorMap.CreateImageTexture(renderer.gl, mat.map.image);
							}
							else {
								colorMap.CreateConstantTexture(renderer.gl, mat.color.toArray());
							}

							let material;
							let shadowMaterials = [];
							let Translation = [transform.modelTransX, transform.modelTransY, transform.modelTransZ];
							let Scale = [transform.modelScaleX, transform.modelScaleY, transform.modelScaleZ];

							let lights = renderer.lights;
							switch (objMaterial) {
								case 'PhongMaterial':
									material = buildPhongMaterial(colorMap, mat.specular.toArray(), lights, "./src/shaders/phongShader/phongVertex.glsl", "./src/shaders/phongShader/phongFragment.glsl");
									for(i=0;i<lights.length;i++){
										let shdowMat = buildShadowMaterial(lights[i].entity, "./src/shaders/shadowShader/shadowVertex.glsl", "./src/shaders/shadowShader/shadowFragment.glsl")
										shadowMaterials.push(
											shdowMat
										);
									}
									break;
							}

							material.then((data) => {
								let meshRender = new MeshRender(renderer.gl, mesh, data);	
								if(autoMove){
									meshRender.setAutoMove();
								}					
								renderer.addMeshRender(meshRender);
							});
							for(let matIndex=0;matIndex<shadowMaterials.length;matIndex++){
								shadowMaterials[matIndex].then((data) => {
									let shadowMeshRender = new MeshRender(renderer.gl, mesh, data);
									renderer.addShadowMeshRender(shadowMeshRender,matIndex);
								});
							}
						}
					});
				}, onProgress, onError);
		});
}
