using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class Cutoff : MonoBehaviour {

	public Transform victim;
	MeshFilter _meshFilter;


	struct VertexInfo{
		public Vector3 vertex;
		public float distanceToPlane;
		public int index;
	}

	void Start()
	{
		_meshFilter = victim.gameObject.GetComponent<MeshFilter> ();
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
				if(vertexInfo[j].distanceToPlane >= 0)
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
				newTriangles.Add (newVertices.Count);
				newVertices.Add (vertexInfo [0].vertex);
				newTriangles.Add (newVertices.Count);
				newVertices.Add (vertexInfo [1].vertex);
				newTriangles.Add (newVertices.Count);
				newVertices.Add (vertexInfo [2].vertex);

			} else if (inPositiveHalfSpaceNum==2) {

				Vector3 newVertexOne = FindContactPointOnEdge(vertexInfo [aNegativeVertexIndex], vertexInfo [(aNegativeVertexIndex + 1)%3], cutPlane, ref edgeDict, ref newGeneratedVertices);
				Vector3 newVertexTwo = FindContactPointOnEdge(vertexInfo [aNegativeVertexIndex], vertexInfo [(aNegativeVertexIndex + 2)%3], cutPlane, ref edgeDict, ref newGeneratedVertices);

				newTriangles.Add (newVertices.Count);
				newVertices.Add (newVertexOne);
				newTriangles.Add (newVertices.Count);
				newVertices.Add (vertexInfo [(aNegativeVertexIndex+1)%3].vertex);
				newTriangles.Add (newVertices.Count);
				newVertices.Add (vertexInfo [(aNegativeVertexIndex+2)%3].vertex);


				newTriangles.Add (newVertices.Count);
				newVertices.Add (newVertexOne);
				newTriangles.Add (newVertices.Count);
				newVertices.Add (vertexInfo [(aNegativeVertexIndex+2)%3].vertex);
				newTriangles.Add (newVertices.Count);
				newVertices.Add (newVertexTwo);


			} else if (inPositiveHalfSpaceNum==1) {
				Vector3 newVertexOne = FindContactPointOnEdge(vertexInfo [aPositiveVertexIndex], vertexInfo [(aPositiveVertexIndex + 1)%3], cutPlane, ref edgeDict, ref newGeneratedVertices);
				Vector3 newVertexTwo = FindContactPointOnEdge(vertexInfo [aPositiveVertexIndex], vertexInfo [(aPositiveVertexIndex + 2)%3], cutPlane, ref edgeDict, ref newGeneratedVertices);


				newTriangles.Add (newVertices.Count);
				newVertices.Add (vertexInfo [aPositiveVertexIndex].vertex);
				newTriangles.Add (newVertices.Count);
				newVertices.Add (newVertexOne);
				newTriangles.Add (newVertices.Count);
				newVertices.Add (newVertexTwo);

			} else if (inPositiveHalfSpaceNum==0) {
				//Whole triangle is culled by plane, just ignore this situation
				continue;
			}
		}

		List<Vector3> contourVertices = FindContour(newGeneratedVertices, cutPlane);
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
		Mesh newMesh = new Mesh ();
		newMesh.vertices = newVertices.ToArray();
		newMesh.triangles = newTriangles.ToArray();
		newMesh.RecalculateNormals ();
		_meshFilter.mesh = newMesh;

	}

	Vector3 FindContactPointOnEdge(VertexInfo vertexa, VertexInfo vertexb, Plane cutPlane, ref Dictionary<string, Vector3> edgeDict, ref List<Vector3> newGeneratedVertices)
	{
		Vector3 newVertex;
		float rayDistance;
		string edgeString;
		if(vertexa.index >= vertexb.index)
			edgeString = vertexa.index.ToString() + vertexb.index.ToString();
		else
			edgeString = vertexb.index.ToString() + vertexa.index.ToString();

		if(edgeDict.ContainsKey(edgeString))
		{
			newVertex = edgeDict[edgeString];
		}
		else
		{
			Ray aRay = new Ray (vertexa.vertex, vertexb.vertex - vertexa.vertex);
			cutPlane.Raycast (aRay, out rayDistance);
			newVertex = aRay.GetPoint (rayDistance);
			newGeneratedVertices.Add (newVertex);
			edgeDict.Add(edgeString, newVertex);
		}
		return newVertex;
	}
	
	List<Vector3> FindContour(List<Vector3> input, Plane plane)
	{
		if(input.Count == 0)
			return new List<Vector3>();
		List<Vector3> output = new List<Vector3> ();
		output.Add (input [0]);
		for(int i = 0; i < input.Count ; i++) {
			Vector3 lastContourPoint = output [output.Count - 1];
			for (int j=0; j<input.Count; j++) 
			{
				if((input[j] - lastContourPoint).magnitude == 0)
					continue;
				if(IsNextContourPoint(lastContourPoint, input[j], input, plane))
				{
					output.Add (input [j]);
					break;
				}
			}

			if((output[0] - output[output.Count - 1]).magnitude == 0)
			{

				return output;
			}
		}


		return output;
	}

	bool IsNextContourPoint(Vector3 lastContourPoint, Vector3 currentPoint, List<Vector3> input, Plane plane)
	{
		for (int i = 0; i<input.Count; i++) 
		{
			if (Vector3.Dot(plane.normal, Vector3.Cross (currentPoint - lastContourPoint, input[i] - lastContourPoint)) < 0) {
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

