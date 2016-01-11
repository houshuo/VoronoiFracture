using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class Cutoff : MonoBehaviour {

	public Transform victim;
	MeshFilter _meshFilter;
	private Mesh _mesh;
	private int[] _triangles;
	private Vector3[] _vertices;


	struct VertexInfo{
		public Vector3 vertex;
		public float distanceToPlane;
		public int index;
	}

	void Initialize()
	{
		_meshFilter = victim.gameObject.GetComponent<MeshFilter> ();
		_mesh = _meshFilter.mesh;
		_triangles = _mesh.triangles;
		_vertices = _mesh.vertices;
	}

	public void Cut()
	{
		Initialize ();
		Plane cutPlane = new Plane (transform.forward, transform.position);
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
		/*List<Vector3> contourVertices = FindContour(newGeneratedVertices, cutPlane);
		if (contourVertices.Count >= 3) {
			for (int i = 1; i<contourVertices.Count - 1; i++) {
				newTriangles.Add (newVertices.Count);
				newVertices.Add (contourVertices[0]);
				newTriangles.Add (newVertices.Count);
				newVertices.Add (contourVertices[i]);
				newTriangles.Add (newVertices.Count);
				newVertices.Add (contourVertices[i+1]);
			}
		}*/
		Mesh newMesh = new Mesh ();
		newMesh.vertices = newVertices.ToArray();
		newMesh.triangles = newTriangles.ToArray();
		_meshFilter.mesh = newMesh;
	}

	Vector3 FindContactPointOnEdge(VertexInfo vertexa, VertexInfo vertexb, Plane cutPlane, ref Dictionary<string, Vector3> edgeDict, ref List<Vector3> newGeneratedVertices)
	{
		Vector3 newVertex;
		float rayDistance;
		string edgeString;
		if(vertexa.index > vertexb.index)
			edgeString = vertexa.index.ToString() + vertexb.index.ToString();
		else
			edgeString = vertexb.index.ToString() + vertexa.index.ToString();

		if(edgeDict.ContainsKey(edgeString))
		{
			newVertex = edgeDict[edgeString];
		}
		else
		{
			Ray firstRay = new Ray (vertexa.vertex, vertexb.vertex - vertexa.vertex);
			cutPlane.Raycast (firstRay, out rayDistance);
			newVertex = firstRay.GetPoint (rayDistance);
			newGeneratedVertices.Add (newVertex);
			edgeDict.Add(edgeString, newVertex);
		}
		return newVertex;
	}
	
	List<Vector3> FindContour(List<Vector3> input, Plane plane)
	{
		List<Vector3> output = new List<Vector3> ();
		output.Add (input [0]);
		input.RemoveAt (0);
		while (input.Count > 0) {
			Vector3 lastContourPoint = output [output.Count - 1];
			for (int i=0; i<input.Count; i++) {
				bool isOnContour = true;
				for (int j = 0; j<input.Count; j++) {
					if (i == j)
						continue;
					if((input[i] - input[j]).magnitude < 0.01)
					{
						input.RemoveAt (j);
						continue;
					}
					if (Vector3.Dot(plane.normal, Vector3.Cross (input[i] - lastContourPoint, input[j] - lastContourPoint)) < 0) {
						isOnContour = false;
						break;
					}
					else
					{
						if( Vector3.Dot (input[i] - lastContourPoint, input[j] - lastContourPoint) > 0)
						{
							if((input[i] - lastContourPoint).magnitude < (input[j] - lastContourPoint).magnitude)
							{
								isOnContour = false;
								break;
							}
							else
							{
								input.RemoveAt (j);
								continue;
							}
						}
					}
				}
				if (isOnContour) {
					output.Add (input [i]);
					input.RemoveAt (i);
					break;
				}
			}
		}
		return output;
	}

	void OnMouseDown()
	{
		Cut ();
	}
}

