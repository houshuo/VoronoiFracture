using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class Cutoff : MonoBehaviour {

	public Transform victim;
	MeshFilter _meshFilter;
	MeshCollider _collider;
	public float epslion = 0.0001f;


	struct VertexInfo{
		public Vector3 vertex;
		public float distanceToPlane;
		public int index;
	}

	void Start()
	{
		_meshFilter = victim.gameObject.GetComponent<MeshFilter> ();
		_collider = victim.gameObject.GetComponent<MeshCollider> ();
	}
	

	public void Cut()
	{
		Mesh _mesh = _meshFilter.mesh;
		int[] _triangles = _mesh.triangles;
		Vector3[] _vertices = _mesh.vertices;
		Vector3 planeNormalInSubspace = victim.InverseTransformDirection (transform.forward);
		Vector3 planePosInSubspace = victim.InverseTransformPoint (transform.position);
		Plane cutPlane = new Plane (planeNormalInSubspace, planePosInSubspace);
		List<Vector3> newVertices = new List<Vector3> ();
		List<Vector3> newGeneratedVertices = new List<Vector3> ();
		List<int> newTriangles = new List<int> ();
		Dictionary<string, Vector3> edgeDict = new Dictionary<string, Vector3> ();
		for (int i = 0; i*3<_triangles.Length; i++) {
			//Get vertices and their distance to plane
			VertexInfo[] vertexInfo = new VertexInfo[3];

			int inPositiveHalfSpaceNum = 0;
			int aPositiveVertexIndex = 0;
			int aNegativeVertexIndex = 0;
			for (int j = 0; j < 3; j++) {
				Vector3 aVertex = _vertices [_triangles[3*i+j]];
				vertexInfo [j] = new VertexInfo();
				vertexInfo[j].distanceToPlane = cutPlane.GetDistanceToPoint (aVertex);
				vertexInfo[j].vertex = aVertex;
				vertexInfo[j].index = _triangles[3*i+j];
				if(vertexInfo[j].distanceToPlane > epslion)
				{
					aPositiveVertexIndex = j;
					inPositiveHalfSpaceNum++;
				}
				else
				{
					aNegativeVertexIndex = j;
				}
			}
			if (inPositiveHalfSpaceNum==3) {
				//Whole triangle is in positive side, just copy the triangle
				AddVertex(vertexInfo [0].vertex, ref newVertices, ref newTriangles);
				AddVertex(vertexInfo [1].vertex, ref newVertices, ref newTriangles);
				AddVertex(vertexInfo [2].vertex, ref newVertices, ref newTriangles);


			} else if (inPositiveHalfSpaceNum==2) {

				Vector3 newVertexOne = FindContactPointOnEdge(vertexInfo [aNegativeVertexIndex], vertexInfo [(aNegativeVertexIndex + 1)%3], cutPlane);
				Vector3 newVertexTwo = FindContactPointOnEdge(vertexInfo [aNegativeVertexIndex], vertexInfo [(aNegativeVertexIndex + 2)%3], cutPlane);

				AddNewGeneratedVertex(newVertexOne, ref newGeneratedVertices);
				AddNewGeneratedVertex(newVertexTwo, ref newGeneratedVertices);

				AddVertex(newVertexOne, ref newVertices, ref newTriangles);
				AddVertex (vertexInfo [(aNegativeVertexIndex+1)%3].vertex, ref newVertices, ref newTriangles);
				AddVertex(vertexInfo [(aNegativeVertexIndex+2)%3].vertex, ref newVertices, ref newTriangles);


				AddVertex(newVertexOne, ref newVertices, ref newTriangles);
				AddVertex(vertexInfo [(aNegativeVertexIndex+2)%3].vertex, ref newVertices, ref newTriangles);
				AddVertex(newVertexTwo, ref newVertices, ref newTriangles);


			} else if (inPositiveHalfSpaceNum==1) {
				Vector3 newVertexOne = FindContactPointOnEdge(vertexInfo [aPositiveVertexIndex], vertexInfo [(aPositiveVertexIndex + 1)%3], cutPlane);
				Vector3 newVertexTwo = FindContactPointOnEdge(vertexInfo [aPositiveVertexIndex], vertexInfo [(aPositiveVertexIndex + 2)%3], cutPlane);
			
				AddNewGeneratedVertex(newVertexOne, ref newGeneratedVertices);
				AddNewGeneratedVertex(newVertexTwo, ref newGeneratedVertices);

				AddVertex(vertexInfo [aPositiveVertexIndex].vertex, ref newVertices, ref newTriangles);
				AddVertex(newVertexOne, ref newVertices, ref newTriangles);
				AddVertex(newVertexTwo, ref newVertices, ref newTriangles);

			} else if (inPositiveHalfSpaceNum==0) {
				//Whole triangle is culled by plane, just ignore this situation
				continue;
			}
		}

		FindContourMesh (newGeneratedVertices, cutPlane, ref newVertices, ref newTriangles);


		Mesh newMesh = new Mesh ();
		newMesh.vertices = newVertices.ToArray();
		newMesh.triangles = newTriangles.ToArray();
		_meshFilter.mesh = newMesh;

		_collider.sharedMesh = newMesh;

	}

	void AddVertex(Vector3 vertex, ref List<Vector3> newVertices, ref List<int> newTriangles)
	{
		int index = newVertices.FindIndex(v => (v - vertex).magnitude < epslion);
		if(index != -1)
		{
			newTriangles.Add (index);
		}
		else
		{
			newTriangles.Add (newVertices.Count);
			newVertices.Add (vertex);
		}
	}

	void AddNewGeneratedVertex(Vector3 vertex, ref List<Vector3> newGeneratedVertices)
	{
		if (newGeneratedVertices.FindIndex (v => (v - vertex).magnitude < epslion) < 0) {
			newGeneratedVertices.Add(vertex);
		}
	}

	Vector3 FindContactPointOnEdge(VertexInfo vertexa, VertexInfo vertexb, Plane cutPlane)
	{
		Vector3 newVertex;
		float rayDistance;
		Ray aRay = new Ray (vertexa.vertex, vertexb.vertex - vertexa.vertex);
		cutPlane.Raycast (aRay, out rayDistance);
		newVertex = aRay.GetPoint (rayDistance);
		return newVertex;
	}
	
	void FindContourMesh(List<Vector3> input, Plane plane, ref List<Vector3> newVertices, ref List<Vector3> newTriangles)
	{
		if(input.Count == 0)
			return new List<Vector3>();
		List<Vector3> contourVertices = new List<Vector3> ();
		contourVertices.Add (input [0]);
		for(int i = 0; i < input.Count ; i++) 
		{
			Vector3 lastContourPoint = contourVertices [contourVertices.Count - 1];
			for (int j=0; j<input.Count; j++) 
			{
				if(input[j] == lastContourPoint)
				{
					continue;
				}
				if(IsNextContourPoint(lastContourPoint, input[j], input, plane))
				{
					contourVertices.Add (input [j]);
					break;
				}
			}
		}
		
		if (contourVertices.Count >= 3) {
			for (int i = 1; i<contourVertices.Count - 1; i++) {
				newTriangles.Add (newVertices.Count);
				newVertices.Add (contourVertices[0]);
				newTriangles.Add (newVertices.Count);
				newVertices.Add (contourVertices[i+1]);
				newTriangles.Add (newVertices.Count);
				newVertices.Add (contourVertices[i]);
			}
		}
		return output;
	}

	bool IsNextContourPoint(Vector3 lastContourPoint, Vector3 currentPoint, List<Vector3> input, Plane plane)
	{
		for (int i = 0; i<input.Count; i++) 
		{
			//Debug.Log (Vector3.Dot(plane.normal, Vector3.Cross (currentPoint - lastContourPoint, input[i] - lastContourPoint)).ToString());
			if (Vector3.Dot(plane.normal, Vector3.Cross (currentPoint - lastContourPoint, input[i] - lastContourPoint)) < -1 * epslion) 
			{
				return false;
			}
		}
		return true;
	}

	void OnMouseDown()
	{
		Cut ();
	}
}

