using UnityEngine;
using System.Linq;
using System.Collections;
using System.Collections.Generic;

public class Cutoff : MonoBehaviour {

	public Transform victim;
	MeshFilter _meshFilter;
	public float epslion = 0.0001f;
	public GameObject newMeshPrefab;


	class VertexInfo{
		public Vector3 vertex;
		public int index;
		public List<int> belongToTriangleIndex;
		public List<string> belongToEdgeIndex;
	}

	class EdgeInfo{
		public int vertexAIndex;
		public int vertexBIndex;
		public List<int> belongToTriangleIndex;
		public int breakVertexIndex;

		public bool IsContainVertex(int vertexIndex)
		{
			return vertexIndex == vertexAIndex || vertexIndex == vertexBIndex;
		}

		public int GetOtherPoint(int vertexIndex)
		{
			if (vertexIndex == vertexAIndex)
				return vertexBIndex;
			else if (vertexIndex == vertexBIndex)
				return vertexAIndex;
			else
				return -1;
		}
	}

	class TriangleInfo{
		public int index;
		public List<string> edges;
		public List<int> vertices;
	}

	void Start()
	{
		_meshFilter = victim.gameObject.GetComponent<MeshFilter> ();
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
		GameObject newMeshLeft = (GameObject)Instantiate (newMeshPrefab, transform.position, transform.rotation);
		MeshFilter newMeshLeftFilter = newMeshLeft.GetComponent<MeshFilter> ();
		List< TriangleInfo> leftTriangles = new List<TriangleInfo>();

		GameObject newMeshRight = (GameObject)Instantiate (newMeshPrefab, transform.position, transform.rotation);
		MeshFilter newMeshRightFilter = newMeshRight.GetComponent<MeshFilter> ();
		List< TriangleInfo> rightTriangles = new List<TriangleInfo>();

		//Get Cutplane
		Vector3 planeNormalInSubspace = victim.InverseTransformDirection (transform.forward);
		Vector3 planePosInSubspace = victim.InverseTransformPoint (transform.position);
		Plane cutPlane = new Plane (planeNormalInSubspace, planePosInSubspace);

		//Get Mesh Info
		Mesh _mesh = _meshFilter.mesh;
		int[] _triangles = _mesh.triangles;
		Vector3[] _vertices = _mesh.vertices;
		Dictionary<int, VertexInfo> vertices = new Dictionary<int, VertexInfo> ();
		Dictionary<string, EdgeInfo> edges = new Dictionary<string, EdgeInfo> ();
		Dictionary<int, TriangleInfo> triangles = new Dictionary<int, TriangleInfo> ();
		for (int i = 0; i<_vertices.Count(); i++) {
			AddVertex(_vertices [i], ref vertices);
		}
		for (int i = 0; i*3<_triangles.Length; i++) {
			VertexInfo[] verticesToBeAdd = new VertexInfo [3];
			for (int j = 0; j < 3; j++) {
				verticesToBeAdd[j] = vertices[_triangles[i*3+j]];
			}
			AddTriangle(verticesToBeAdd, ref edges, ref triangles);
		}

		//Cut Mesh
		int triangleNum = triangles.Count;
		for(int i = 0; i < triangleNum; i++)
		{
			int inPositiveHalfSpaceNum = 0;
			int aPositiveVertexIndex = 0;
			int aNegativeVertexIndex = 0;
			TriangleInfo triangle = triangles[i];
			for(int j = 0; j < 3; j++)
			{
				VertexInfo avertex = vertices[triangle.vertices[j]];
				float distanceToPlane = cutPlane.GetDistanceToPoint (avertex.vertex);
				if(distanceToPlane > 0)
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
				leftTriangles.Add (triangle);
				
			} else if (inPositiveHalfSpaceNum==2) {
				List<EdgeInfo> crossEdge = new List<EdgeInfo>();
				foreach(string edgeString in triangle.edges)
				{
					EdgeInfo edge = edges[edgeString];
					if(edge.IsContainVertex(triangle.vertices[aNegativeVertexIndex]))
					{
						crossEdge.Add (edge);
					}
				}

				for(int k = 0; k < 2; k++)
				{
					EdgeInfo edge = crossEdge[k];
					if(edge.breakVertexIndex == -1)
					{
						Vector3 breakPoint = FindContactPointOnEdge(vertices[edge.vertexAIndex], vertices[edge.vertexBIndex], cutPlane);
						VertexInfo vertex = AddVertex(breakPoint, ref vertices);
						edge.breakVertexIndex = vertex.index;
					}
				}

				VertexInfo[] triangleA = new VertexInfo[3]{vertices[crossEdge[0].breakVertexIndex], 
					vertices[crossEdge[0].GetOtherPoint(triangle.vertices[aNegativeVertexIndex])], 
					vertices[crossEdge[1].GetOtherPoint(triangle.vertices[aNegativeVertexIndex])]};

				VertexInfo[] triangleB = new VertexInfo[3]{vertices[crossEdge[0].breakVertexIndex], 
					vertices[crossEdge[1].GetOtherPoint(triangle.vertices[aNegativeVertexIndex])], 
					vertices[crossEdge[1].breakVertexIndex]};

				VertexInfo[] triangleC = new VertexInfo[3]{vertices[triangle.vertices[aNegativeVertexIndex]], 
					vertices[crossEdge[0].breakVertexIndex],
					vertices[crossEdge[1].breakVertexIndex]};

				leftTriangles.Add (AddTriangle(triangleA, ref edges, ref triangles));
				leftTriangles.Add (AddTriangle(triangleB, ref edges, ref triangles));
				rightTriangles.Add (AddTriangle(triangleC, ref edges, ref triangles));
				
			} else if (inPositiveHalfSpaceNum==1) {
				List<EdgeInfo> crossEdge = new List<EdgeInfo>();
				foreach(string edgeString in triangle.edges)
				{
					EdgeInfo edge = edges[edgeString];
					if(edge.IsContainVertex(triangle.vertices[aPositiveVertexIndex]))
					{
						crossEdge.Add (edge);
					}
				}
				
				for(int k = 0; k < 2; k++)
				{
					EdgeInfo edge = crossEdge[k];

					if(edge.breakVertexIndex == -1)
					{
						Vector3 breakPoint= FindContactPointOnEdge(vertices[edge.vertexAIndex], vertices[edge.vertexBIndex], cutPlane);
						VertexInfo vertex = AddVertex(breakPoint, ref vertices);
						edge.breakVertexIndex = vertex.index;
					}
				}

				VertexInfo[] triangleA = new VertexInfo[3]{vertices[crossEdge[0].breakVertexIndex], 
					vertices[crossEdge[1].GetOtherPoint(triangle.vertices[aPositiveVertexIndex])], 
					vertices[crossEdge[0].GetOtherPoint(triangle.vertices[aPositiveVertexIndex])]};
				
				VertexInfo[] triangleB = new VertexInfo[3]{vertices[crossEdge[0].breakVertexIndex], 
					vertices[crossEdge[1].breakVertexIndex],
					vertices[crossEdge[1].GetOtherPoint(triangle.vertices[aPositiveVertexIndex])]};

				
				VertexInfo[] triangleC = new VertexInfo[3]{vertices[triangle.vertices[aPositiveVertexIndex]], 
					vertices[crossEdge[0].breakVertexIndex],
					vertices[crossEdge[1].breakVertexIndex]};
				
				rightTriangles.Add (AddTriangle(triangleA, ref edges, ref triangles));
				rightTriangles.Add (AddTriangle(triangleB, ref edges, ref triangles));
				leftTriangles.Add (AddTriangle(triangleC, ref edges, ref triangles));

			} else if (inPositiveHalfSpaceNum==0) {
				//Whole triangle is culled by plane, just ignore this situation
				rightTriangles.Add (triangle);
			}
		}
		//Cut mesh end
		List<Vector3> leftVertices = new List<Vector3> ();
		Dictionary<int, int> verticeIndexCorrespondingDict = new Dictionary<int, int> ();
		List<int> leftTriangleIndex = new List<int> ();
		for (int i = 0; i < leftTriangles.Count; i++) {
			TriangleInfo triangle = leftTriangles[i];
			foreach(int vIndex in triangle.vertices)
			{
				VertexInfo vertexInfo = vertices[vIndex];
				if(!verticeIndexCorrespondingDict.ContainsKey(vertexInfo.index))
				{
					verticeIndexCorrespondingDict.Add(vertexInfo.index, leftVertices.Count);
					leftVertices.Add (vertexInfo.vertex);

				}
				leftTriangleIndex.Add (verticeIndexCorrespondingDict[vertexInfo.index]);
			}
		}

		newMeshLeftFilter.mesh.vertices = leftVertices.ToArray();
		newMeshLeftFilter.mesh.triangles = leftTriangleIndex.ToArray();
		newMeshLeftFilter.mesh.RecalculateNormals ();

		verticeIndexCorrespondingDict.Clear ();
		List<Vector3> rightVertices = new List<Vector3> ();
		List<int> rightTriangleIndex = new List<int> ();
		for (int i = 0; i < rightTriangles.Count; i++) {
			TriangleInfo triangle = rightTriangles[i];
			foreach(int vIndex in triangle.vertices)
			{
				VertexInfo vertexInfo = vertices[vIndex];
				if(!verticeIndexCorrespondingDict.ContainsKey(vertexInfo.index))
				{
					verticeIndexCorrespondingDict.Add(vertexInfo.index, rightVertices.Count);
					rightVertices.Add (vertexInfo.vertex);
					
				}
				rightTriangleIndex.Add (verticeIndexCorrespondingDict[vertexInfo.index]);
			}
		}

		newMeshRightFilter.mesh.vertices = rightVertices.ToArray();
		newMeshRightFilter.mesh.triangles = rightTriangleIndex.ToArray();
		newMeshRightFilter.mesh.RecalculateNormals ();

		Destroy (victim.gameObject);
	}

	VertexInfo AddVertex(Vector3 pos, ref Dictionary<int, VertexInfo> vertices)
	{
		VertexInfo vertex = new VertexInfo();
		vertex.vertex = pos;
		vertex.index = vertices.Count;
		vertex.belongToTriangleIndex = new List<int> ();
		vertex.belongToEdgeIndex = new List<string> ();
		vertices.Add (vertex.index, vertex);
		return vertex;
	}

	TriangleInfo AddTriangle(VertexInfo[] verticesToAdd, ref Dictionary<string, EdgeInfo> edges, ref Dictionary<int, TriangleInfo> triangles)
	{
		TriangleInfo triangle = new TriangleInfo();
		triangle.vertices = new List<int> ();
		triangle.edges = new List<string> ();
		triangle.index = triangles.Count;
		for (int i = 0; i<3; i++) {
			verticesToAdd[i].belongToTriangleIndex.Add (triangle.index);
			triangle.vertices.Add(verticesToAdd[i].index);
			EdgeInfo edge;
			string edgeString = GetEdgeString (verticesToAdd[i].index, verticesToAdd[(i+1)%3].index);
			if (edges.ContainsKey (edgeString))
				edge = edges [edgeString];
			else {
				edge = new EdgeInfo ();
				edge.belongToTriangleIndex = new List<int> ();
				edge.vertexAIndex = verticesToAdd[i].index;
				edge.vertexBIndex = verticesToAdd[(i+1)%3].index;
				edge.breakVertexIndex = -1;
				edges.Add (edgeString, edge);
			}
			triangle.edges.Add(edgeString);
			edge.belongToTriangleIndex.Add (triangle.index);
			verticesToAdd[i].belongToEdgeIndex.Add(edgeString);
		}
		triangles.Add (triangle.index, triangle);
		return triangle;
	}


	#region FindContactPoint
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
	#endregion

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

