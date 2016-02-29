//Purpose: A file that holds the code that students fill in


function areaAroundPoint(p, verticies){

  var total = 0;
  for (var i = 0; i < verticies.length; i++){
    var a = verticies[i];
    if(i == verticies.length - 1){
      var b = verticies[0];
    } else {
      var b = verticies[i + 1];
    }
    var pa = vec3.create();
    var pb = vec3.create();
    var cross = vec3.create();

    vec3.subtract(pa, a, p);
    vec3.subtract(pb, b, p);
    vec3.cross(cross, pa, pb);
    total += vec3.length(cross) / 2;
    }
    return total;
 }
function projVector(u, v) {
    var scale = (vec3.dot(u, v)/vec3.dot(v, v));//The scale in front of v is (u dot v) / (v dot v)
    var projv = vec3.create(); //Allocate a vector to hold the output
    vec3.scale(projv, v, scale); //Scale v by the appropriate amount
    return projv; //Return the result
}
function inPolygon(p, verticies){

  var pArea = areaAroundPoint(p, verticies);
  var tArea = areaAroundPoint(verticies[0], verticies);

  if (Math.abs(pArea-tArea) <= 1e-4){
    return true;
  }
  return false;
}
function arraysAreIdentical(arr1, arr2){
    if (arr1.length !== arr2.length) return false;
    for (var i = 0, len = arr1.length; i < len; i++){
        if (arr1[i] !== arr2[i]){
            return false;
        }
    }
    return true;
}
//Given a ray described by an initial point P0 and a direction V both in
//world coordinates, check to see
//if it intersects the polygon described by "vertices," an array of vec3
//values describing the location of the polygon vertices in its child frame.
//mvMatrix is a matrix describing how to transform "vertices" into world coordinates
//which you will have to do to get the correct intersection in world coordinates.
//Be sure to compute the plane normal only after you have transformed the points,
//and be sure to only compute intersections which are inside of the polygon
//(you can assume that all polygons are convex and use the area method)
function rayIntersectPolygon(P0, V, vertices, mvMatrix) {
    //TODO: Fill this in
    //Step 1: Make a new array of vec3s which holds "vertices" transformed
    //to world coordinates (hint: vec3 has a function "transformMat4" which is useful)

    var tVerticies = [];
    for (var i = 0; i < vertices.length;i++){
      var tVertex = vec3.create();
      vec3.transformMat4(tVertex, vertices[i], mvMatrix);
      tVerticies.push(tVertex);
    }
    var bP0 = vec3.create();
    vec3.subtract(bP0, P0, tVerticies[1]);
    //Step 2: Compute the plane normal of the plane spanned by the transformed vertices
    var ba = vec3.create();
    var bc = vec3.create();
    var norm = vec3.create();

    vec3.subtract(ba, tVerticies[0], tVerticies[1]);
    vec3.subtract(bc, tVerticies[2], tVerticies[1]);

    vec3.cross(norm, ba, bc);
    //Step 3: Perform ray intersect plane
    var denom = vec3.dot(V, norm);
    var numer = -1 * vec3.dot(bP0, norm);

    if(denom == 0){
      return null;
    }

    var tIntersect = numer / denom;

    if(tIntersect <= 0){
        return null;
    }

    var P = vec3.create();
    vec3.scaleAndAdd(P, P0, V, tIntersect);

    //Step 4: Check to see if the intersection point is inside of the transformed polygon
    //You can assume that the polygon is convex.  If you use the area test, you can
    //allow for some wiggle room in the two areas you're comparing (e.g. absolute difference
    //not exceeding 1e-4)

    if(!inPolygon(P, tVerticies)){
        return null;
    }

    //Step 5: Return the intersection point if it exists or null if it's outside
    //of the polygon or if the ray is perpendicular to the plane normal (no intersection)

    return {t:tIntersect, P:P};
}
function caluculateT(endNode, startNode, V){
    var endNodeT = vec3.create();
    vec3.subtract(endNodeT,endNode.pos,startNode.pos);

    if(V[0] != 0){
        var scalar = (1/V[0]);
        var t = endNodeT[0]*scalar;
    } else if(V[1] != 0){
        var scalar = (1/V[1]);
        var t = endNodeT[1]*scalar;
    } else{
        var scalar = (1/V[2]);
        var t = endNodeT[2]*scalar;
    }
    return t;
}
function convertToWorldCoordinate(vertices,mvMatrix){
    var rVertices = [];
    for (var i = 0; i < vertices.length; i++) {
        var rVec = vec3.create();
        vec3.transformMat4(rVec, vertices[i], mvMatrix)
        rVertices.push(rVec);
    }
    return rVertices;
}
//Returns the box that surrounds the given mesh object in world coordinates
function boxMesh(mesh,mvMatrix){
    var minX = Infinity;
    var minY = Infinity;
    var minZ = Infinity;
    var maxX = Number.NEGATIVE_INFINITY;
    var maxY = Number.NEGATIVE_INFINITY;
    var maxZ = Number.NEGATIVE_INFINITY;

    for (var f = 0; f < mesh.faces.length; f++) {
        var worldVertices = convertToWorldCoordinate(mesh.faces[f].getVerticesPos(),mvMatrix);
        for(var i = 0; i < worldVertices.length;i++){
            if(worldVertices[i][0] > maxX){ maxX = worldVertices[i][0];}
            if(worldVertices[i][1] > maxY){ maxY = worldVertices[i][1];}
            if(worldVertices[i][2] > maxZ){ maxZ = worldVertices[i][2];}

            if(worldVertices[i][0] < minX){ minX = worldVertices[i][0];}
            if(worldVertices[i][1] < minY){ minY = worldVertices[i][1];}
            if(worldVertices[i][2] < minZ){ minZ = worldVertices[i][2];}
        }
    }
    return {minX:minX, minY:minY, minZ:minZ, maxX:maxX, maxY:maxY, maxZ:maxZ};
}
function mergeBoxes(box1,box2){
    var minX = 0;
    var minY = 0;
    var minZ = 0;
    var maxX = 0;
    var maxY = 0;
    var maxZ = 0;

    if(box1.maxX > box2.maxX){ maxX = box1.maxX;} else {maxX = box2.maxX;}
    if(box1.maxY > box2.maxY){ maxY = box1.maxY;} else {maxY = box2.maxY;}
    if(box1.maxZ > box2.maxZ){ maxZ = box1.maxZ;} else {maxZ = box2.maxZ;}

    if(box1.minX < box2.minX){ minX = box1.minX;} else {minX = box2.minX;}
    if(box1.minY < box2.minY){ minY = box1.minY;} else {minY = box2.minY;}
    if(box1.minZ < box2.minZ){ minZ = box1.minZ;} else {minZ = box2.minZ;}
    return {minX:minX, minY:minY, minZ:minZ, maxX:maxX, maxY:maxY, maxZ:maxZ};
}

/// Global Variables ///
var p = 0;
function addImageSourcesFunctions(scene) {
    scene.preComputeBoxes = function(node,mvMatrix){
        var returnBox;
        if(!("children" in node)){
            //compute bounding box for this node
            var box = boxMesh(node.mesh,mvMatrix);
            //Add dummy node as box and child is the node
            returnBox = {
              childrenNodes: [node],
              childrenBoxes: null,
              box: box,
            };
            return returnBox;
        }
        else{
            var returnBox = {
              childrenNodes: [],
              childrenBoxes: [],
              box: {minX:Infinity, minY:Infinity, minZ:Infinity, maxX:Number.NEGATIVE_INFINITY, maxY:Number.NEGATIVE_INFINITY, maxZ:Number.NEGATIVE_INFINITY},
            };
            //Compute subboxes and merge with master box
            for(var i = 0; i < node.children.length;i++){

                //Find new mvMatrix
                var nextmvMatrix = mat4.create();
                mat4.mul(nextmvMatrix, mvMatrix, node.children[i].transform);

                //Find next childBox
                var childBox = scene.preComputeBoxes(node.children[i],nextmvMatrix);

                //Update current box
                returnBox.box = mergeBoxes(returnBox.box,childBox.box);
                returnBox.childrenBoxes.push(childBox);
            }
        }
        return returnBox;
    }
    //Setup all of the functions that students fill in that operate directly
    //on the scene

    //Purpose: A recursive function provided which helps to compute intersections of rays
    //with all faces in the scene, taking into consideration the scene graph structure
    //Inputs: P0 (vec3): Ray starting point, V (vec3): ray direction
    //node (object): node in scene tree to process,
    //mvMatrix (mat4): Matrix to put geometry in this node into world coordinates
    //excludeFace: Pointer to face object to be excluded (don't intersect with
    //the face that this point lies on)
    //Returns: null if no intersection,
    //{tmin:minimum t along ray, PMin(vec3): corresponding point, faceMin:Pointer to mesh face hit first}

    //NOTE: Calling this function with node = scene and an identity matrix for mvMatrix
    //will start the recursion at the top of the scene tree in world coordinates
    scene.rayIntersectFaces = function(P0, V, node, mvMatrix, excludeFace) {
        var tmin = Infinity;//The parameter along the ray of the nearest intersection
        var PMin = null;    //The point of intersection corresponding to the nearest interesection
        var faceMin = null; //The face object corresponding to the nearest intersection

        if (node === null) {
            return null;
        }
        if ('mesh' in node) { //Make sure it's not just a dummy transformation node
            var mesh = node.mesh;
            for (var f = 0; f < mesh.faces.length; f++) {
                if (mesh.faces[f] == excludeFace) {
                    continue;//Don't re-intersect with the face this point lies on
                }
                //Intersect the ray with this polygon
                var res = rayIntersectPolygon(P0, V, mesh.faces[f].getVerticesPos(), mvMatrix);
                if (!(res === null) && (res.t < tmin)) {
                    tmin = res.t;
                    PMin = res.P;
                    faceMin = mesh.faces[f];
                }
            }
        }

        if ('children' in node) {
            //Recursively check the meshes of the children to make sure the ray
            //doesn't intersect any of them first
            for (var i = 0; i < node.children.length; i++) {
                var nextmvMatrix = mat4.create();
                //Multiply on the right by the next transformation of the child
                //node
                mat4.mul(nextmvMatrix, mvMatrix, node.children[i].transform);
                //Recursively intersect with the child node
                var cres = scene.rayIntersectFaces(P0, V, node.children[i], nextmvMatrix, excludeFace);
                if (!(cres === null) && (cres.tmin < tmin)) {
                    tmin = cres.tmin;
                    PMin = cres.PMin;
                    faceMin = cres.faceMin;
                }
            }
        }
        if (PMin === null) {
            return null;
        }
        return {tmin:tmin, PMin:PMin, faceMin:faceMin};
    }

    scene.computeImageSourcesHelper = function(node,order,obj,mvMatrix){
        if('mesh' in node){
            for(var m = 0; m < node.mesh.faces.length; m++){
                //Current objs genFace != current face
                if(obj.genFace == null || !(obj.genFace == node.mesh.faces[m])){
                    var norm = vec3.create();
                    var ba = vec3.create();
                    var bc = vec3.create();
                    var pa = vec3.create();
                    var reflect_pt = vec3.create(); // Return pt
                    var vertices = convertToWorldCoordinate(node.mesh.faces[m].getVerticesPos(),mvMatrix);

                    vec3.subtract(ba,vertices[0],vertices[1]);
                    vec3.subtract(bc,vertices[2],vertices[1]);
                    vec3.subtract(pa,vertices[0], obj.pos); // P is the source in world coordinates
                    vec3.cross(norm,ba,bc);

                    var proj = projVector(pa,norm);
                    vec3.scaleAndAdd(reflect_pt,obj.pos,proj,2);

                    scene.imsources.push({
                      pos: reflect_pt,
                      order: order,
                      rcoeff: node.rcoeff,
                      parent: obj,
                      genFace: node.mesh.faces[m],
                      mvMatrix: mvMatrix
                  });
              }
            }
        }
        if('children' in node){
            for (var i = 0; i < node.children.length; i++) {
                var nextmvMatrix = mat4.create();
                mat4.mul(nextmvMatrix, mvMatrix, node.children[i].transform);

                scene.computeImageSourcesHelper(node.children[i],order,obj,nextmvMatrix);
            }

        }
    }
    //Purpose: Fill in the array scene.imsources[] with a bunch of source
    //objects.  It's up to you what you put in the source objects, but at
    //the very least each object needs a field "pos" describing its position
    //in world coordinates so that the renderer knows where to draw it
    //You will certainly also need to save along pointers from an image source
    //to its parent so that when you trace paths back you know where to aim
    //Recursion is highly recommended here, since you'll be making images of
    //images of images (etc...) reflecting across polygon faces.

    //Inputs: order (int) : The maximum number of bounces to take
    scene.computeImageSources = function(order) {
        scene.source.order = 0;//Store an order field to figure out how many
        //bounces a particular image represents
        scene.source.rcoeff = 1.0;//Keep track of the reflection coefficient of the node that
        //gave rise to this source
        scene.source.parent = null;//Keep track of the image source's parent
        scene.source.genFace = null;//Keep track of the mesh face that generated this image
        scene.source.mvMatrix = null;
        //Remember not to reflect an image across the face that just generated it,
        //or you'll get its parent image.  This information can also be used later
        //when tracing back paths
        scene.imsources = [scene.source];

        //TODO: Fill the rest of this in.  Be sure to reflect images across faces
        //in world coordinates, not the faces in the original mesh coordinates
        //See the "rayIntersectFaces" function above for an example of how to loop
        //through faces in a mesh
        console.log(scene);
        var mvMatrix = [1,0,0,0 ,0,1,0,0 ,0,0,1,0 ,0,0,0,1];

        for(var p = 1; p <= order;p++){
            for(var i = scene.imsources.length-1; i >= 0; i--) {
                if(scene.imsources[i].order == p-1){
                    scene.computeImageSourcesHelper(scene,p,scene.imsources[i],mvMatrix);
                }
            }
        }
        console.log("# of sources:",scene.imsources.length);
    }

    //Purpose: Based on the extracted image sources, trace back paths from the
    //receiver to the source, checking to make sure there are no occlusions
    //along the way.  Remember, you're always starting by tracing a path from
    //the receiver to the image, and then from the intersection point with
    //that image's corresponding face to the image's parent, and so on
    //all the way until you get back to the original source.

    //Fill in the array scene.paths, where each element of the array is itself
    //an array of objects describing vertices along the path, starting
    //with the receiver and ending with the source.  Each object in each path
    //array should contain a field "pos" which describes the position, as well
    //as an element "rcoeff" which stores the reflection coefficient at that
    //part of the path, which will be used to compute decays in "computeInpulseResponse()"
    //Don't forget the direct path from source to receiver!
    scene.extractPathsHelper = function(startNode,endNode,subpath,impulse,totalDist) {
        var mvMatrix = [1,0,0,0 ,0,1,0,0 ,0,0,1,0 ,0,0,0,1];

        var V = vec3.create();
        vec3.subtract(V,startNode.pos,endNode.pos);

        //Break condition -- no bouncing necessary since it connects with the source
        if(startNode == scene.source){
            var directIntersect = scene.rayIntersectFaces(endNode.pos, V, scene, mvMatrix, endNode.genFace);
            //var t = caluculateT(endNode,startNode,V);

            // Check if there is a direct intersection between startNode and endNode excluding the generation face of the endNode
            if(directIntersect == null || directIntersect.tmin >= 1){
                //Add subpath
                subpath.push(scene.source);
                scene.paths.push(subpath);

                //Add impulse material
                var d = vec3.distance(scene.source.pos,endNode.pos);
                imp = impulse *  (1/(1+Math.pow(d,p)));
                scene.impulses.push([imp,d+totalDist]);
            }
            return;
        }
        //Calculate t from endNode to startNode generation face
        var tToBouncePt = rayIntersectPolygon(endNode.pos, V, startNode.genFace.getVerticesPos(), startNode.mvMatrix);
        //Calculate minimum t from endNode to startNode
        var intersect = scene.rayIntersectFaces(endNode.pos, V, scene, mvMatrix, startNode.genFace);

        if((tToBouncePt != null && tToBouncePt.t > 0 && tToBouncePt.t <= 1) && (intersect == null || intersect.tmin > tToBouncePt.t)){ // t must be greater than 0, since its a ray

                var intermediateNode = {
                  pos: tToBouncePt.P,
                  order: startNode.order,
                  rcoeff: startNode.rcoeff,
                  parent: startNode.parent,
                  genFace: startNode.genFace,
                  mvMatrix: startNode.mvMatrix
                };
                subpath.push(intermediateNode);

                var d = vec3.distance(startNode.pos,endNode.pos);
                var imp = impulse *  1/(1+Math.pow(d,p)) * startNode.rcoeff;

                return scene.extractPathsHelper(intermediateNode.parent,intermediateNode,subpath,imp,totalDist+d);
        }
        else{ return; }
        return;
    }
    scene.extractPaths = function() {

        scene.paths = [];
        //Identity Matrix
        var mvMatrix = [1,0,0,0 ,0,1,0,0 ,0,0,1,0 ,0,0,0,1];

        //compute boxes
        var boxes = scene.preComputeBoxes(scene,mvMatrix);
        console.log("boxes");
        console.log(boxes);
        //First check the direct path from the source
        var V = vec3.create();
        var path = [scene.receiver];
        vec3.subtract(V,scene.receiver.pos,scene.source.pos);
        var intersect = scene.rayIntersectFaces(scene.source.pos, V, scene, mvMatrix, null);

        //If intersect = null there are no intersections
        if(intersect == null || intersect.tmin > 1){
            path.push(scene.source);
            scene.paths.push(path);

            var d = vec3.distance(scene.receiver.pos,scene.source.pos);
            var imp = 1/(1+Math.pow(d,p));
            scene.impulses.push([imp,d]);
        }
        // Recursively check bounces for other nodes
        for(var q = 1; q < scene.imsources.length; q++){
            scene.extractPathsHelper(scene.imsources[q],scene.receiver,[scene.receiver],1,0);
        }

        //TODO: Finish this. Extract the rest of the paths by backtracing from
        //the image sources you calculated.  Return an array of arrays in
        //scene.paths.  Recursion is highly recommended
        //Each path should start at the receiver and end at the source
        //(or vice versa), so scene.receiver should be the first element
        //and scene.source should be the last element of every array in
        //scene.paths
    }


    //Inputs: Fs: Sampling rate (samples per second)
    scene.computeImpulseResponse = function(Fs) {
        var SVel = 340;//Sound travels at 340 meters/second
        var Fs = 44100;
        //TODO: Finish this.  Be sure to scale each bounce by 1/(1+r^p),
        //where r is the length of the line segment of that bounce in meters
        //and p is some integer less than 1 (make it smaller if you want the
        //paths to attenuate less and to be more echo-y as they propagate)
        //Also be sure to scale by the reflection coefficient of each material
        //bounce (you should have stored this in extractPaths() if you followed
        //those directions).  Use some form of interpolation to spread an impulse
        //which doesn't fall directly in a bin to nearby bins
        scene.impulsesResp = [];
        //Save the result into the array
        //Impulse Response Start/end
        for(var j = 0; j < scene.impulses.length;j++){
            //Direct Path
            //convert to time
            var sampleNum = Math.floor(Fs * scene.impulses[j][1] / SVel);
            while(sampleNum >= scene.impulseResp.length){
                scene.impulseResp.push(0);
            }
            scene.impulseResp[sampleNum] += scene.impulses[j][0];
        }
    }
}
