//Purpose: A file that holds the code that students fill in


function areaAroundPoint(p, verticies){

  var total = 0;
  for (i = 0; i < verticies.length; i++){
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

  if (abs(pArea-tArea) <= 1e-4){
    return true;
  }
  return false;
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
    for (vertex in vertices){
      var tVertex = vec3.create();
      vec3.transformMat4(tVertex, vertex, mvMatrix);
    }

    //Step 2: Compute the plane normal of the plane spanned by the transformed vertices
    var ba = vec3.create();
    var bc = vec3.create();
    var norm = vec3.create();

    vec3.subtract(ba, tVerticies[0], tVerticies[1]);
    vec3.subtract(bc, tVerticies[2], tVerticies[1]);

    vec3.cross(norm, bc, ba);

    //Step 3: Perform ray intersect plane
    var denom;
    var numer;

    denom = vec3.dot(V, norm);
    numer = -1 * vec3.dot(P0, norm);

    if(demon == 0){
      return null;
    }

    var tIntersect = numer / denom;
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

    return {t:tIntersect, P:P}; //These are dummy values, but you should return
    //both an intersection point and a parameter t.  The parameter t will be used to sort
    //intersections in order of occurrence to figure out which one happened first
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
function createReflections(order,obj,source){
    var sources = [];
    console.log(obj);
    for(var m = 0; m < obj.mesh.faces.length;m++){                         // Reflect over everyface in the
        if(obj.parent == null || obj.genFace != obj.mesh.faces[m]){ // don't go back over the same face!!
            if(order == 2){
                console.log("faces");
                console.log(obj.genFace);
                console.log(obj.mesh.faces[m]);
                console.log("facesEnd");
            }
            var norm = vec3.create();
            var ba = vec3.create();
            var bc = vec3.create();
            var pa = vec3.create();
            var reflect_pt = vec3.create(); // Return pt
            var vertices = convertToWorldCoordinate(obj.mesh.faces[m].getVerticesPos(),obj.transform);
            vec3.subtract(ba,vertices[0],vertices[1]);
            vec3.subtract(bc,vertices[2],vertices[1]);
            vec3.subtract(pa,vertices[0], source.pos); // P is the source in world coordinates
            vec3.cross(norm,ba,bc);
            var proj = projVector(pa,norm);
            vec3.scaleAndAdd(reflect_pt,source.pos,proj,2);
            sources.push({
              pos: reflect_pt,
              order: order,
              rcoeff: obj.mesh.rcoeff,
              parent: obj,
              genFace: obj.mesh.faces[m]
          });
      }
    }
    if ('children' in obj) {
        for(var x = 0; x < obj.children.length;x++){
            if(obj.children[x].order == order-1){
                sources.push.apply(sources,createReflections(order,obj.children[x],source));
            }
        }
    }
    return sources;
}
function addImageSourcesFunctions(scene) {
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
        var PMin = null;//The point of intersection corresponding to the nearest interesection
        var faceMin = null;//The face object corresponding to the nearest intersection
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
        //Remember not to reflect an image across the face that just generated it,
        //or you'll get its parent image.  This information can also be used later
        //when tracing back paths
        scene.imsources = [scene.source];

        //TODO: Fill the rest of this in.  Be sure to reflect images across faces
        //in world coordinates, not the faces in the original mesh coordinates
        //See the "rayIntersectFaces" function above for an example of how to loop
        //through faces in a mesh

        console.log(scene.source);

        for(var p = 1; p <= order;p++){
            for(var i = 0; i < scene.children.length; i++) {
                if(!('order' in scene.children[i]) || scene.children[i].order == p-1){

                    scene.imsources.push.apply(scene.imsources,createReflections(p,scene.children[i],scene.source));
                }
            }
        }
        console.log(scene.imsources);
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
    scene.extractPaths = function() {
        scene.paths = [];

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
        //TODO: Finish this.  Be sure to scale each bounce by 1/(1+r^p),
        //where r is the length of the line segment of that bounce in meters
        //and p is some integer less than 1 (make it smaller if you want the
        //paths to attenuate less and to be more echo-y as they propagate)
        //Also be sure to scale by the reflection coefficient of each material
        //bounce (you should have stored this in extractPaths() if you followed
        //those directions).  Use some form of interpolation to spread an impulse
        //which doesn't fall directly in a bin to nearby bins
        //Save the result into the array scene.impulseResp[]
    }
}
