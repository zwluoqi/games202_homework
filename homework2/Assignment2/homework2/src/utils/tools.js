const { MathUtils } = require("../../lib/three");

function getRotationPrecomputeL(precompute_L, rotationMatrix){
	let newprecompute_L = Array.from(precompute_L);

	var M33 = computeSquareMatrix_3by3(rotationMatrix);
	let sh2 = [];
	sh2.push(precompute_L[1]);
	sh2.push(precompute_L[2]);
	sh2.push(precompute_L[3]);
	var newsh2 = math.multiply(M33,sh2);
	newprecompute_L[1] = newsh2._data[0];
	newprecompute_L[2] = newsh2._data[1];
	newprecompute_L[3] = newsh2._data[2];

	var M55 = computeSquareMatrix_5by5(rotationMatrix);
	var sh3 = [];
	sh3.push(precompute_L[4]);
	sh3.push(precompute_L[5]);
	sh3.push(precompute_L[6]);
	sh3.push(precompute_L[7]);
	sh3.push(precompute_L[8]);
	var newsh3 = math.multiply(M55,sh3);
	newprecompute_L[4] = newsh3._data[0];
	newprecompute_L[5] = newsh3._data[1];
	newprecompute_L[6] = newsh3._data[2];
	newprecompute_L[7] = newsh3._data[3];
	newprecompute_L[8] = newsh3._data[4];

	return newprecompute_L;
}

function computeSquareMatrix_3by3(rotationMatrix){ // 计算方阵SA(-1) 3*3 
	
	// 1、pick ni - {ni}
	let n1 = [1, 0, 0, 0]; let n2 = [0, 0, 1, 0]; let n3 = [0, 1, 0, 0];

	// 2、{P(ni)} - A  A_inverse
	var p1 = SHEval(n1[0],n1[1],n1[2],3);
	var p2 = SHEval(n2[0],n2[1],n2[2],3);
	var p3 = SHEval(n3[0],n3[1],n3[2],3);
	var A = math.matrix([[p1[1],p1[2],p1[3]],
		[p2[1],p2[2],p2[3]],
		[p3[1],p3[2],p3[3]]]
		);
	var A_Inverse = math.inv(A);


	// 3、用 R 旋转 ni - {R(ni)}
	var r1 = math.multiply(rotationMatrix,  n1);
	var r2 = math.multiply(rotationMatrix,  n2);
	var r3 = math.multiply(rotationMatrix,  n3);

	// 4、R(ni) SH投影 - S

	r1 = r1._data;
	r2 = r2._data;
	r3 = r3._data;

	var pr1 = SHEval(r1[0],r1[1],r1[2],3);
	var pr2 = SHEval(r2[0],r2[1],r2[2],3);
	var pr3 = SHEval(r3[0],r3[1],r3[2],3);
	var S = math.matrix([[pr1[1],pr1[2],pr1[3]],
		[pr2[1],pr2[2],pr2[3]],
		[pr3[1],pr3[2],pr3[3]]]
		);



	// 5、S*A_inverse
	var M = math.multiply(S,A_Inverse);
	return M;
}

function computeSquareMatrix_5by5(rotationMatrix){ // 计算方阵SA(-1) 5*5
	
	// 1、pick ni - {ni}
	let k = 1 / math.sqrt(2);
	let n1 = [1, 0, 0, 0]; let n2 = [0, 0, 1, 0]; let n3 = [k, k, 0, 0]; 
	let n4 = [k, 0, k, 0]; let n5 = [0, k, k, 0];

	// 2、{P(ni)} - A  A_inverse
	var p1 = SHEval(n1[0],n1[1],n1[2],3);
	var p2 = SHEval(n2[0],n2[1],n2[2],3);
	var p3 = SHEval(n3[0],n3[1],n3[2],3);
	var p4 = SHEval(n4[0],n4[1],n4[2],3);
	var p5 = SHEval(n5[0],n5[1],n5[2],3);

	var A = math.matrix([
		[p1[4],p1[5],p1[6],p1[7],p1[8]],
		[p2[4],p2[5],p2[6],p2[7],p2[8]],
		[p3[4],p3[5],p3[6],p3[7],p3[8]],
		[p4[4],p4[5],p4[6],p4[7],p4[8]],
		[p5[4],p5[5],p5[6],p5[7],p5[8]],
	]
		);
	var A_Inverse = math.inv(A);


	// 3、用 R 旋转 ni - {R(ni)}
	var r1 = math.multiply(rotationMatrix,  n1);
	var r2 = math.multiply(rotationMatrix,  n2);
	var r3 = math.multiply(rotationMatrix,  n3);
	var r4 = math.multiply(rotationMatrix,  n4);
	var r5 = math.multiply(rotationMatrix,  n5);

	// 4、R(ni) SH投影 - S
	r1 = r1._data;
	r2 = r2._data;
	r3 = r3._data;
	r4 = r4._data;
	r5 = r5._data;

	var pr1 = SHEval(r1[0],r1[1],r1[2],3);
	var pr2 = SHEval(r2[0],r2[1],r2[2],3);
	var pr3 = SHEval(r3[0],r3[1],r3[2],3);
	var pr4 = SHEval(r4[0],r4[1],r4[2],3);
	var pr5 = SHEval(r5[0],r5[1],r5[2],3);

	var S =  math.matrix([
		[pr1[4],pr1[5],pr1[6],pr1[7],pr1[8]],
		[pr2[4],pr2[5],pr2[6],pr2[7],pr2[8]],
		[pr3[4],pr3[5],pr3[6],pr3[7],pr3[8]],
		[pr4[4],pr4[5],pr4[6],pr4[7],pr4[8]],
		[pr5[4],pr5[5],pr5[6],pr5[7],pr5[8]],
	]
		);


	// 5、S*A_inverse
	var M = math.multiply(S,A_Inverse);
	return M;
}

function mat4Matrix2mathMatrix(rotationMatrix){

	let mathMatrix = [];
	for(let i = 0; i < 4; i++){
		let r = [];
		for(let j = 0; j < 4; j++){
			r.push(rotationMatrix[i*4+j]);
		}
		mathMatrix.push(r);
	}
	return math.matrix(mathMatrix)

}

function getMat3ValueFromRGB(precomputeL){

    let colorMat3 = [];
    for(var i = 0; i<3; i++){
        colorMat3[i] = mat3.fromValues( precomputeL[0][i], precomputeL[1][i], precomputeL[2][i],
										precomputeL[3][i], precomputeL[4][i], precomputeL[5][i],
										precomputeL[6][i], precomputeL[7][i], precomputeL[8][i] ); 
	}
    return colorMat3;
}