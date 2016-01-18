using UnityEngine;
using System.Linq;
using System.Collections;
using System.Collections.Generic;

public class Cutoff : MonoBehaviour {

	public Transform victim;
	MeshFilter _meshFilter;
	public float epslion = 0.0001f;
	public GameObject newMeshPrefab;

	void OnMouseDown()
	{
		Cut ();
	}


	class VertexInfo{
		public Vector3 vertex;
		public Vector3 normal;
		public int index;
		public List<int> belongToTriangleIndex;
		public List<string> belongToEdgeIndex;
	}

	class EdgeInfo{
		public int vertexAIndex;
		public int vertexBIndex;
		public List<int> belongToTriangleIndex;
		public bool isContourEdge;
		public bool IsCrossEdge(EdgeInfo edge)
		{
			return true;
		}

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

		public int GetOtherTriangle(int triangleIndex)
		{
			if (belongToTriangleIndex.Count < 2)
				return -1;
			else if (triangleIndex == belongToTriangleIndex [0])
				return belongToTriangleIndex [1];
			else if (triangleIndex == belongToTriangleIndex [1])
				return belongToTriangleIndex [0];
			else
				return -1;
		}

		public string GetSelfEdgeString()
		{
			return GetEdgeString (vertexAIndex, vertexBIndex);
		}

		public static string GetEdgeString(int vertexAIndex, int vertexBIndex)
		{
			if (vertexAIndex > vertexBIndex)
				return vertexAIndex.ToString () + "-" + vertexBIndex.ToString ();
			else if (vertexBIndex > vertexAIndex)
				return vertexBIndex.ToString () + "-" + vertexAIndex.ToString ();
			else
				return "";
		}
	}

	class TriangleInfo{
		public int index;
		public List<string> edges;
		public List<int> vertices;
	}

	class DelaunayTriangleInfo: TriangleInfo
	{
		public static bool IsAPointInsideTrianglesCircumcircle(Vector3 point, Vector3 vertexA, Vector3 vertexB, Vector3 vertexC)
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
	}

	void Start()
	{
		_meshFilter = victim.gameObject.GetComponent<MeshFilter> ();
	}



	public void Cut()
	{
		//Get Mesh
		GameObject newMeshLeft = (GameObject)Instantiate (newMeshPrefab, victim.position, victim.rotation);
		MeshFilter newMeshLeftFilter = newMeshLeft.GetComponent<MeshFilter> ();
		List< TriangleInfo> leftTriangles = new List<TriangleInfo>();

		GameObject newMeshRight = (GameObject)Instantiate (newMeshPrefab, victim.position, victim.rotation);
		MeshFilter newMeshRightFilter = newMeshRight.GetComponent<MeshFilter> ();
		List< TriangleInfo> rightTriangles = new List<TriangleInfo>();

		//Get Cutplane
		Vector3 planeNormalInSubspace = victim.InverseTransformDirection (transform.forward);
		Vector3 planePosInSubspace = victim.InverseTransformPoint (transform.position);
		Plane cutPlane = new Plane (planeNormalInSubspace, planePosInSubspace);

		//New Generated Points
		Dictionary<string, VertexInfo> newGeneratedPoints = new Dictionary<string, VertexInfo> ();

		//Get Mesh Info
		Mesh _mesh = _meshFilter.mesh;
		int[] _triangles = _mesh.triangles;
		Vector3[] _vertices = _mesh.vertices;
		Vector3[] _normals = _mesh.normals;
		Dictionary<int, VertexInfo> vertices = new Dictionary<int, VertexInfo> ();
		Dictionary<string, EdgeInfo> edges = new Dictionary<string, EdgeInfo> ();
		Dictionary<int, TriangleInfo> triangles = new Dictionary<int, TriangleInfo> ();
		for (int i = 0; i<_vertices.Count(); i++) {
			AddVertex(_vertices [i], _normals[i], ref vertices);
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

				for(int k = 0; k < 2; k++)
				{
					string edgeString = EdgeInfo.GetEdgeString(triangle.vertices[aNegativeVertexIndex], triangle.vertices[(aNegativeVertexIndex+1+k)%3]);
					EdgeInfo edge = edges[edgeString];
					if(!newGeneratedPoints.ContainsKey(edgeString))
					{
						Vector3 breakPoint;
						Vector3 normal;
						FindContactPointAndNormalOnEdge(vertices[edge.vertexAIndex], vertices[edge.vertexBIndex], cutPlane, out breakPoint, out normal);
						VertexInfo vertex = AddVertex(breakPoint, normal, ref vertices);
						newGeneratedPoints[edgeString] = vertex;
					}
					crossEdge.Add (edge);
				}

				VertexInfo[] triangleA = new VertexInfo[3]{newGeneratedPoints[crossEdge[0].GetSelfEdgeString()], 
					vertices[crossEdge[0].GetOtherPoint(triangle.vertices[aNegativeVertexIndex])], 
					vertices[crossEdge[1].GetOtherPoint(triangle.vertices[aNegativeVertexIndex])]};

				VertexInfo[] triangleB = new VertexInfo[3]{newGeneratedPoints[crossEdge[0].GetSelfEdgeString()], 
					vertices[crossEdge[1].GetOtherPoint(triangle.vertices[aNegativeVertexIndex])], 
					newGeneratedPoints[crossEdge[1].GetSelfEdgeString()]};

				VertexInfo[] triangleC = new VertexInfo[3]{vertices[triangle.vertices[aNegativeVertexIndex]], 
					newGeneratedPoints[crossEdge[0].GetSelfEdgeString()],
					newGeneratedPoints[crossEdge[1].GetSelfEdgeString()]};

				leftTriangles.Add (AddTriangle(triangleA, ref edges, ref triangles));
				leftTriangles.Add (AddTriangle(triangleB, ref edges, ref triangles));
				rightTriangles.Add (AddTriangle(triangleC, ref edges, ref triangles));
				
			} else if (inPositiveHalfSpaceNum==1) {
				List<EdgeInfo> crossEdge = new List<EdgeInfo>();
				for(int k = 0; k < 2; k++)
				{
					string edgeString = EdgeInfo.GetEdgeString(triangle.vertices[aPositiveVertexIndex], triangle.vertices[(aPositiveVertexIndex+1+k)%3]);
					EdgeInfo edge = edges[edgeString];

					if(!newGeneratedPoints.ContainsKey(edgeString))
					{
						Vector3 breakPoint;
						Vector3 normal;
						FindContactPointAndNormalOnEdge(vertices[edge.vertexAIndex], vertices[edge.vertexBIndex], cutPlane, out breakPoint, out normal);
						VertexInfo vertex = AddVertex(breakPoint, normal, ref vertices);
						newGeneratedPoints[edgeString] = vertex;
					}
					crossEdge.Add (edge);
				}

				VertexInfo[] triangleA = new VertexInfo[3]{newGeneratedPoints[crossEdge[0].GetSelfEdgeString()], 
					vertices[crossEdge[0].GetOtherPoint(triangle.vertices[aPositiveVertexIndex])], 
					vertices[crossEdge[1].GetOtherPoint(triangle.vertices[aPositiveVertexIndex])]};
				
				VertexInfo[] triangleB = new VertexInfo[3]{newGeneratedPoints[crossEdge[0].GetSelfEdgeString()], 
					vertices[crossEdge[1].GetOtherPoint(triangle.vertices[aPositiveVertexIndex])],
					newGeneratedPoints[crossEdge[1].GetSelfEdgeString()]};

				
				VertexInfo[] triangleC = new VertexInfo[3]{vertices[triangle.vertices[aPositiveVertexIndex]], 
					newGeneratedPoints[crossEdge[0].GetSelfEdgeString()],
					newGeneratedPoints[crossEdge[1].GetSelfEdgeString()]};
				
				rightTriangles.Add (AddTriangle(triangleA, ref edges, ref triangles));
				rightTriangles.Add (AddTriangle(triangleB, ref edges, ref triangles));
				leftTriangles.Add (AddTriangle(triangleC, ref edges, ref triangles));

			} else if (inPositiveHalfSpaceNum==0) {
				//Whole triangle is culled by plane, just ignore this situation
				rightTriangles.Add (triangle);
			}
		}
		//Cut mesh end

		//Add new face
		List<VertexInfo> newGenerateFaceContour = FindContourVertex (newGeneratedPoints.Values.ToList(), cutPlane.normal);
		List<VertexInfo> newGenerateFaceLeft = new List<VertexInfo> ();
		List<VertexInfo> newGenerateFaceRight = new List<VertexInfo> ();
		foreach (VertexInfo avertex in newGenerateFaceContour) {
			newGenerateFaceLeft.Add (AddVertex (avertex.vertex, cutPlane.normal, ref vertices));
			newGenerateFaceRight.Add (AddVertex (avertex.vertex, -1*cutPlane.normal, ref vertices));
		}
		if (newGenerateFaceContour.Count > 2) {
			for(int i = 1; i<newGenerateFaceContour.Count -1; i++)
			{
				VertexInfo[] triangleLeft = new VertexInfo[3]{newGenerateFaceLeft[0],
					newGenerateFaceLeft[i],
					newGenerateFaceLeft[i+1]};

				VertexInfo[] triangleRight = new VertexInfo[3]{newGenerateFaceRight[0],
					newGenerateFaceRight[i+1],
					newGenerateFaceRight[i]};

				leftTriangles.Add (AddTriangle(triangleLeft, ref edges, ref triangles));
				rightTriangles.Add (AddTriangle(triangleRight, ref edges, ref triangles));
			}
		}

		//Add new face end

		//re-assemble mesh
		List<Vector3> leftVertices = new List<Vector3> ();
		List<Vector3> leftNormals = new List<Vector3> ();
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
					leftNormals.Add (vertexInfo.normal);
				}
				leftTriangleIndex.Add (verticeIndexCorrespondingDict[vertexInfo.index]);
			}
		}

		newMeshLeftFilter.mesh.vertices = leftVertices.ToArray();
		newMeshLeftFilter.mesh.normals = leftNormals.ToArray ();
		newMeshLeftFilter.mesh.triangles = leftTriangleIndex.ToArray();
		newMeshLeftFilter.mesh.RecalculateNormals ();
		newMeshLeftFilter.mesh.RecalculateBounds ();

		verticeIndexCorrespondingDict.Clear ();
		List<Vector3> rightVertices = new List<Vector3> ();
		List<Vector3> rightNormals = new List<Vector3> ();
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
					rightNormals.Add (vertexInfo.normal);
				}
				rightTriangleIndex.Add (verticeIndexCorrespondingDict[vertexInfo.index]);
			}
		}

		newMeshRightFilter.mesh.vertices = rightVertices.ToArray();
		newMeshRightFilter.mesh.normals = rightNormals.ToArray ();
		newMeshRightFilter.mesh.triangles = rightTriangleIndex.ToArray();
		newMeshRightFilter.mesh.RecalculateNormals ();
		newMeshRightFilter.mesh.RecalculateBounds ();

		Destroy (victim.gameObject);
	}

	VertexInfo AddVertex(Vector3 pos, Vector3 normal, ref Dictionary<int, VertexInfo> vertices)
	{
		VertexInfo vertex = new VertexInfo();
		vertex.vertex = pos;
		vertex.normal = normal;
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
			string edgeString = EdgeInfo.GetEdgeString (verticesToAdd[i].index, verticesToAdd[(i+1)%3].index);
			if (edges.ContainsKey (edgeString))
				edge = edges [edgeString];
			else {
				edge = new EdgeInfo ();
				edge.belongToTriangleIndex = new List<int> ();
				edge.vertexAIndex = verticesToAdd[i].index;
				edge.vertexBIndex = verticesToAdd[(i+1)%3].index;
				edges.Add (edgeString, edge);
			}
			triangle.edges.Add(edgeString);
			edge.belongToTriangleIndex.Add (triangle.index);
			verticesToAdd[i].belongToEdgeIndex.Add(edgeString);
		}
		triangles.Add (triangle.index, triangle);
		return triangle;
	}

	EdgeInfo AddEdge(VertexInfo vertexa, VertexInfo vertexb, ref Dictionary<string, EdgeInfo> edges)
	{
		string edgeString = EdgeInfo.GetEdgeString(vertexa.index, vertexb.index);
		EdgeInfo edge;
		if (!edges.ContainsKey (edgeString)) {
			edge = new EdgeInfo ();
			edge.vertexAIndex = vertexa.index;
			edge.vertexBIndex = vertexb.index;
			edges [edgeString] = edge;
		} else {
			edge = edges[edgeString];
		}
		return edge;
	}

	void FindContactPointAndNormalOnEdge(VertexInfo vertexa, VertexInfo vertexb, Plane cutPlane, out Vector3 pos, out Vector3 normal)
	{
		float rayDistance;
		Ray aRay = new Ray (vertexa.vertex, vertexb.vertex - vertexa.vertex);
		cutPlane.Raycast (aRay, out rayDistance);
		pos = aRay.GetPoint (rayDistance);

		float ratio = rayDistance / (vertexb.vertex - vertexa.vertex).magnitude;
		normal = Vector3.Lerp (vertexa.vertex, vertexb.vertex, ratio);
	}

	#region FindContourMesh
	List<VertexInfo> FindContourVertex(List<VertexInfo> contour, Vector3 planeNormal)
	{
		if (contour.Count == 0)
			return contour;
		contour = contour.OrderBy(v => v.vertex.x).ThenBy(v => v.vertex.y).ToList();
		List<VertexInfo> duplicated = new List<VertexInfo> ();
		for (int i = 0; i < contour.Count; i++) {
			if((contour[i].vertex - contour[(i+1)%(contour.Count-1)].vertex).magnitude < epslion)
				duplicated.Add (contour[i]);
		}

		contour.RemoveAll (v => duplicated.Contains (v));

		for(int i = 0; i < contour.Count - 1; i++) 
		{
			VertexInfo lastContourPointIndex = contour [i];
			for (int j=i+1; j<contour.Count; j++) 
			{
				if(IsNextContourPoint(lastContourPointIndex, contour[j], contour, planeNormal))
				{
					VertexInfo nextContourPoint = contour[j];
					contour.RemoveAt(j);
					contour.Insert(i+1, nextContourPoint);
					break;
				}
			}
		}
		return contour;

	}

	bool IsNextContourPoint(VertexInfo lastContourPoint, VertexInfo currentPoint, List<VertexInfo> contourPointSet, Vector3 planeNormal)
	{
		for (int i = 0; i<contourPointSet.Count; i++) 
		{
			//Debug.Log (Vector3.Dot(plane.normal, Vector3.Cross (currentPoint - lastContourPoint, input[i] - lastContourPoint)).ToString());
			if (Vector3.Dot(planeNormal, Vector3.Cross (currentPoint.vertex - lastContourPoint.vertex, contourPointSet[i].vertex - lastContourPoint.vertex)) >  epslion) 
			{
				return false;
			}
		}
		return true;
	}
	#endregion
	
	#region Delaunay
	void FindContourMeshDelaunay(List<VertexInfo> contour, Plane plane)
	{
		contour = contour.OrderBy(v=>v.vertex.x).ThenBy(v => v.vertex.z).ToList();
		List<VertexInfo> duplicated = new List<VertexInfo> ();
		for (int i = 0; i < contour.Count; i++) {
			if((contour[i].vertex - contour[(i+1)%(contour.Count-1)].vertex).magnitude < epslion)
				duplicated.Add (contour[i]);
		}
		contour.RemoveAll (v => duplicated.Contains (v));

		Dictionary<string, EdgeInfo> edges = new Dictionary<string, EdgeInfo> ();

		DelaunayDivideAndConquer (contour, ref edges);
	}
	
	void DelaunayDivideAndConquer(List<VertexInfo> vertices, ref Dictionary<string, EdgeInfo> edges)
	{
		if (vertices.Count == 2) {
			EdgeInfo edge = AddEdge (vertices [0], vertices [1], ref edges);
			edge.isContourEdge = true;
		} else if (vertices.Count == 3) {
			for (int i = 0; i < 3; i++) {
				EdgeInfo edge = AddEdge (vertices [i], vertices [(i + 1) % 2], ref edges);
				edge.isContourEdge = true;
			}
		} else if (vertices.Count > 3) {
			Dictionary<string, EdgeInfo> leftEdges = new Dictionary<string, EdgeInfo>();
			Dictionary<string, EdgeInfo> rightEdges = new Dictionary<string, EdgeInfo>();
			DelaunayDivideAndConquer(vertices.Take(vertices.Count/2).ToList(), ref leftEdges);
			DelaunayDivideAndConquer(vertices.Skip(vertices.Count/2).ToList(), ref rightEdges);
		}
	}
	#endregion
}

