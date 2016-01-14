using UnityEngine;
using System.Linq;
using System.Collections;
using System.Collections.Generic;

public class Cutoff : MonoBehaviour {

	public Transform victim;
	MeshFilter _meshFilter;
	MeshCollider _collider;
	public float epslion = 0.0001f;
	public GameObject newMeshPrefab;


	struct VertexInfo{
		public Vector3 vertex;
		public int index;
		public List<int> belongToTriangleIndex;
		public List<int> belongToEdgeIndex;
	}

	struct EdgeInfo{
		public int vertexAIndex;
		public int vertexBIndex;
		public List<int> belongToTriangleIndex;
	}

	struct TriangleInfo{
		public int index;
		public List<int> Vertices;
	}

	void Start()
	{
		_meshFilter = victim.gameObject.GetComponent<MeshFilter> ();
		_collider = victim.gameObject.GetComponent<MeshCollider> ();
	}

	private string GetEdgeString(int vertexAIndex, int vertexBIndex)
	{
		if (vertexAIndex > vertexBIndex)
			return vertexAIndex.ToString () + "-" + vertexBIndex.ToString ();
		else if (vertexBIndex > vertexAIndex)
			return vertexBIndex.ToString () + "-" + vertexAIndex.ToString ();
		else
			return "";
	}

	public void Cut()
	{
		//Get Mesh
		Mesh _mesh = _meshFilter.mesh;
		int[] _triangles = _mesh.triangles;
		Vector3[] _vertices = _mesh.vertices;

		GameObject newMeshLeft = Instantiate (newMeshPrefab, transform.position, transform.rotation);
		GameObject newMeshRight = Instantiate (newMeshPrefab, transform.position, transform.rotation);

		//Get Cutplane
		Vector3 planeNormalInSubspace = victim.InverseTransformDirection (transform.forward);
		Vector3 planePosInSubspace = victim.InverseTransformPoint (transform.position);
		Plane cutPlane = new Plane (planeNormalInSubspace, planePosInSubspace);

		//Get Mesh Info
		Dictionary<int, VertexInfo> vertices = new Dictionary<int, VertexInfo> ();
		Dictionary<string, EdgeInfo> edges = new Dictionary<string, EdgeInfo> ();
		Dictionary<int, TriangleInfo> triangles = new Dictionary<int, TriangleInfo> ();
		for (int i = 0; i*3<_triangles.Length; i++) {
			TriangleInfo triangle = new TriangleInfo ();
			triangle.index = triangles.Count;
			triangle.Vertices = new List<int> ();
			triangles.Add (triangle.index, triangle);
			for (int j = 0; j < 3; j++) {
				VertexInfo vertex;
				if (vertices.ContainsKey (_triangles [3 * i + j])) {
					vertex = vertices [_triangles [3 * i + j]];
				} else {
					vertex = new VertexInfo ();
					vertex.vertex = _vertices [_triangles [3 * i + j]];
					vertex.index = _triangles [3 * i + j];
					vertex.belongToTriangleIndex = new List<int> ();
					vertex.belongToEdgeIndex = new List<int> ();
					vertices.Add (vertex.index, vertex);
				}
				vertex.belongToTriangleIndex.Add (triangle.index);
				triangle.Vertices.Add (vertex.index);

				int nextVertexIndex = _triangles [3 * i + (j + 1) % 3];
				string edgeString = GetEdgeString (vertex.index, nextVertexIndex);
				EdgeInfo edge;
				if (edges.Contains (edgeString))
					edge = edges [edgeString];
				else {
					edge = new EdgeInfo ();
					edge.belongToTriangleIndex = new List<int> ();
					edge.vertexAIndex = _triangles [3 * i + j];
					edge.vertexBIndex = _triangles [3 * i + (j + 1) % 3];
					edges.Add (edgeString, edge);
				}
				edge.belongToTriangleIndex.Add (triangle.index);
			}
		}
		for(int i = 0; i < triangles.Count; i++)
		{
			int inPositiveHalfSpaceNum = 0;
			int aPositiveVertexIndex = 0;
			int aNegativeVertexIndex = 0;
			TriangleInfo triangle = triangles[i];
			for(int j = 0; j < 3; j++)
			{
				VertexInfo avertex = vertices[triangle.Vertices[j]];
				float distanceToPlane = cutPlane.GetDistanceToPoint (avertex.vertex);
				if(distanceToPlane > 0)
				{
					aPositiveVertexIndex = avertex.index;
					inPositiveHalfSpaceNum++;
				}
				else
				{
					aNegativeVertexIndex = avertex.index
				}
			}

			if (inPositiveHalfSpaceNum==3) {
				//Whole triangle is in positive side, just copy the triangle, Do nothing
				
			} else if (inPositiveHalfSpaceNum==2) {
				
				Vector3 newVertexOne = FindContactPointOnEdge(vertexInfo [aNegativeVertexIndex], vertexInfo [(aNegativeVertexIndex + 1)%3], cutPlane);
				Vector3 newVertexTwo = FindContactPointOnEdge(vertexInfo [aNegativeVertexIndex], vertexInfo [(aNegativeVertexIndex + 2)%3], cutPlane);
				
				AddVertex(newVertexOne, ref newVertices, ref newTriangles);
				AddVertex (vertexInfo [(aNegativeVertexIndex+1)%3].vertex, ref newVertices, ref newTriangles);
				AddVertex(vertexInfo [(aNegativeVertexIndex+2)%3].vertex, ref newVertices, ref newTriangles);
				
				
				AddVertex(newVertexOne, ref newVertices, ref newTriangles);
				AddVertex(vertexInfo [(aNegativeVertexIndex+2)%3].vertex, ref newVertices, ref newTriangles);
				AddVertex(newVertexTwo, ref newVertices, ref newTriangles);
				
				AddNewGeneratedVertex(newVertexOne, ref newGeneratedVertices);
				AddNewGeneratedVertex(newVertexTwo, ref newGeneratedVertices);
				
				
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
		int index = newVertices.FindIndex(v => v == vertex);
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
	
	void FindContourMesh(List<Vector3> input, Plane plane, ref List<Vector3> newVertices, ref List<int> newTriangles)
	{
		if (input.Count == 0)
			return;
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

	#region Delaunay
	bool IsAPointInsideTrianglesCircumcircle(Vector3 point, Vector3 vertexA, Vector3 vertexB, Vector3 vertexC)
	{
		Vector3 edgeAB = vertexA - vertexB;
		Vector3 edgeBC = vertexB - vertexC;
		Vector3 edgeCA = vertexC - vertexA;
		float denominator = 2 * (Vector3.Cross(edgeAB, edgeBC)).magnitude;
		float radius = edgeAB.magnitude * edgeBC.magnitude * edgeCA.magnitude / denominator;

		float alpha = Vector3.Dot (edgeBC, edgeBC) * Vector3.Dot(edgeAB , -1 * edgeCA) / denominator;
		float beta = Vector3.Dot (edgeCA, edgeCA) * Vector3.Dot(-1 * edgeAB, edgeBC) / denominator;
		float gamma = Vector3.Dot(edgeAB, edgeAB) * Vector3.Dot(edgeCA, -1 * edgeBC) / denominator;

		Vector3 center = vertexA * alpha + vertexB * beta + vertexC * gamma;

		return (point - center).magnitude < radius;
	}

	/*void FindContourMeshDelaunay(List<Vector3> input, Plane plane, ref List<Vector3> newVertices, ref List<int> newTriangles)
	{
		if (input.Count == 0)
			return;
		IEnumerable<Vector3> sortedVerticesEnumerator = input.OrderBy(v=>v.x).ThenBy(v => v.z);
		input = sortedVerticesEnumerator.ToList ();
		
		DelaunayDivideAndConquer (input, vertices, ref newVertices, ref );
	}
	
	void DelaunayDivideAndConquer(List<int> input, List<Vector3> vertices, ref List<Vector3> newVertices, ref List<int> newTriangles)
	{
		if (input.Count == 3) {
			
		}
	}*/
	#endregion
}

