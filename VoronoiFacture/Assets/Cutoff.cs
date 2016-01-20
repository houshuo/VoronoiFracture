using UnityEngine;
using System.Linq;
using System;
using System.Collections;
using System.Collections.Generic;

public class Cutoff : MonoBehaviour {

	public Transform victim;
	MeshFilter _meshFilter;
	public float epslion = 0.0001f;
	public GameObject newMeshPrefab;

	private Plane cutPlane;
	private Dictionary<int, VertexInfo> vertices;
	private Dictionary<string, EdgeInfo> edges;
	private Dictionary<int, TriangleInfo> triangles;
	private List<VertexInfo> contourVertices;

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
		//initialize public variables
		vertices = new Dictionary<int, VertexInfo> ();
		edges = new Dictionary<string, EdgeInfo> ();
		triangles = new Dictionary<int, TriangleInfo> ();
		//Get Cutplane
		Vector3 planeNormalInSubspace = victim.InverseTransformDirection (transform.forward);
		Vector3 planePosInSubspace = victim.InverseTransformPoint (transform.position);
		cutPlane = new Plane (planeNormalInSubspace, planePosInSubspace);

		//Get Mesh
		GameObject newMeshLeft = (GameObject)Instantiate (newMeshPrefab, victim.position, victim.rotation);
		MeshFilter newMeshLeftFilter = newMeshLeft.GetComponent<MeshFilter> ();
		List< TriangleInfo> leftTriangles = new List<TriangleInfo>();

		GameObject newMeshRight = (GameObject)Instantiate (newMeshPrefab, victim.position, victim.rotation);
		MeshFilter newMeshRightFilter = newMeshRight.GetComponent<MeshFilter> ();
		List< TriangleInfo> rightTriangles = new List<TriangleInfo>();

		//New Generated Points
		Dictionary<string, VertexInfo> newGeneratedPoints = new Dictionary<string, VertexInfo> ();

		//Get Mesh Info
		Mesh _mesh = _meshFilter.mesh;
		int[] _triangles = _mesh.triangles;
		Vector3[] _vertices = _mesh.vertices;
		Vector3[] _normals = _mesh.normals;

		for (int i = 0; i<_vertices.Count(); i++) {
			AddVertex(_vertices [i], _normals[i]);
		}
		for (int i = 0; i*3<_triangles.Length; i++) {
			VertexInfo[] verticesToBeAdd = new VertexInfo [3];
			for (int j = 0; j < 3; j++) {
				verticesToBeAdd[j] = vertices[_triangles[i*3+j]];
			}
			AddTriangle(verticesToBeAdd);
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
						VertexInfo vertex = AddVertex(breakPoint, normal);
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

				leftTriangles.Add (AddTriangle(triangleA));
				leftTriangles.Add (AddTriangle(triangleB));
				rightTriangles.Add (AddTriangle(triangleC));
				
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
						VertexInfo vertex = AddVertex(breakPoint, normal);
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
				
				rightTriangles.Add (AddTriangle(triangleA));
				rightTriangles.Add (AddTriangle(triangleB));
				leftTriangles.Add (AddTriangle(triangleC));

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
			newGenerateFaceLeft.Add (AddVertex (avertex.vertex, cutPlane.normal));
			newGenerateFaceRight.Add (AddVertex (avertex.vertex, -1*cutPlane.normal));
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

				leftTriangles.Add (AddTriangle(triangleLeft));
				rightTriangles.Add (AddTriangle(triangleRight));
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

	VertexInfo AddVertex(Vector3 pos, Vector3 normal)
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

	TriangleInfo AddTriangle(VertexInfo[] verticesToAdd)
	{
		if (!IsRightOfVector (verticesToAdd [0], verticesToAdd [1], verticesToAdd [2])) {
			VertexInfo tmp = verticesToAdd[2];
			verticesToAdd[2] = verticesToAdd[1];
			verticesToAdd[1] = tmp;
		}
		TriangleInfo triangle = new TriangleInfo();
		triangle.vertices = new List<int> ();
		triangle.edges = new List<string> ();
		triangle.index = triangles.Count;
		for (int i = 0; i<3; i++) {
			verticesToAdd[i].belongToTriangleIndex.Add (triangle.index);
			triangle.vertices.Add(verticesToAdd[i].index);
			EdgeInfo edge = _AddEdge (verticesToAdd[i], verticesToAdd[i+1]);
			triangle.edges.Add(edge.GetSelfEdgeString());
			edge.belongToTriangleIndex.Add (triangle.index);
			verticesToAdd[i].belongToEdgeIndex.Add(edge.GetSelfEdgeString());
		}
		triangles.Add (triangle.index, triangle);
		return triangle;
	}

	EdgeInfo _AddEdge(VertexInfo vertexa, VertexInfo vertexb)
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
		vertexa.belongToEdgeIndex.Add (edgeString);
		vertexb.belongToEdgeIndex.Add (edgeString);
		return edge;
	}

	EdgeInfo AddEdge(VertexInfo vertexa, VertexInfo vertexb)
	{
		List<int> vertexaAdjancentVertices = vertexa.belongToEdgeIndex.Select (eIndex => edges [eIndex].GetOtherPoint (vertexa.index)).ToList();
		List<int> vertexbAdjancentVertices = vertexb.belongToEdgeIndex.Select (eIndex => edges [eIndex].GetOtherPoint (vertexb.index)).ToList();
		List<int> commonVertex = vertexaAdjancentVertices.Intersect (vertexbAdjancentVertices).ToList ();
		if (commonVertex.Count > 0) {
			foreach(int vIndex in commonVertex)
			{
				VertexInfo[] newTriangle = new VertexInfo[3]{vertexa, vertexb, vertices[vIndex]};
			}
			return edges[EdgeInfo.GetEdgeString(vertexa.index, vertexb.index)];
		} else {
			return _AddEdge(vertexa, vertexb);
		}
	}

	string RemoveEdge(VertexInfo vertexa, VertexInfo vertexb)
	{
		string edgeString = EdgeInfo.GetEdgeString(vertexa.index, vertexb.index);
		if (edges.ContainsKey (edgeString)) {
			EdgeInfo edge = edges[edgeString];
			foreach(int triangle in edge.belongToTriangleIndex)
			{
				RemoveTriangle(triangles[triangle]);
			}

			vertexa.belongToEdgeIndex.Remove(edgeString);
			vertexb.belongToEdgeIndex.Remove (edgeString);

			edges.Remove (edgeString);
		}
		return edgeString;
	}

	void RemoveTriangle(TriangleInfo triangle)
	{
		foreach (int vIndex in triangle.vertices) {
			vertices[vIndex].belongToTriangleIndex.Remove(triangle.index);
		}
		foreach (string eIndex in triangle.edges) {
			edges[eIndex].belongToTriangleIndex.Remove(triangle.index);
		}
		triangles.Remove (triangle.index);
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
	}
	
	Dictionary<string, EdgeInfo>  DelaunayDivideAndConquer(List<VertexInfo> vertices)
	{
		Dictionary<string, EdgeInfo> newGeneratedEdges = new Dictionary<string, EdgeInfo> ();
		if (vertices.Count == 2) {
			EdgeInfo edge = AddEdge (vertices [0], vertices [1]);
			newGeneratedEdges.Add (edge.GetSelfEdgeString(), edge);
		} else if (vertices.Count == 3) {
			for (int i = 0; i < 3; i++) {
				EdgeInfo edge = AddEdge (vertices [i], vertices [(i + 1) % 2]);
				newGeneratedEdges.Add (edge.GetSelfEdgeString(), edge);
			}
		} else if (vertices.Count > 3) {
			List<VertexInfo> leftVertices = vertices.Take(vertices.Count/2).ToList();
			List<VertexInfo> rightVertices = vertices.Skip(vertices.Count/2).ToList();
			Dictionary<string, EdgeInfo> leftEdges = DelaunayDivideAndConquer(leftVertices);
			Dictionary<string, EdgeInfo> rightEdges = DelaunayDivideAndConquer(rightVertices);
			KeyValuePair<VertexInfo, VertexInfo> lowBoundEdge = FindHullEdge(leftVertices, leftEdges, rightVertices, rightEdges, true);
			KeyValuePair<VertexInfo, VertexInfo> upperBoundEdge = FindHullEdge(leftVertices, leftEdges, rightVertices, rightEdges, false);
			VertexInfo L = lowBoundEdge.Key;
			VertexInfo R = lowBoundEdge.Value;
			while(lowBoundEdge.Key != upperBoundEdge.Key && lowBoundEdge.Value != upperBoundEdge.Value)
			{
				bool A = false;
				bool B = false;
				EdgeInfo edge = AddEdge(L, R);
				newGeneratedEdges.Add (edge.GetSelfEdgeString(), edge);
				VertexInfo R1 = FindPrevVertex(R, L, false);
				if(IsRightOfVector(R, L, R1))
				{
					VertexInfo R2 = FindPrevVertex(R, R1, false);
					while(!TriangleInfo.IsAPointInsideTrianglesCircumcircle(R2.vertex, R1.vertex, L.vertex, R.vertex))
					{
						string edgeString = RemoveEdge(R, R1);
						newGeneratedEdges.Remove(edgeString);
						VertexInfo tmp = R2;
						R2 = FindPrevVertex(R, R1, false);
						R1 = tmp;
					}
				}
				else
				{
					A = true;
				}

				VertexInfo L1 = FindPrevVertex(L, R, true);
				if(IsRightOfVector(R, L, L1))
				{
					VertexInfo L2 = FindPrevVertex(L, L1, true);
					while(!TriangleInfo.IsAPointInsideTrianglesCircumcircle(L2.vertex, L1.vertex, L.vertex, R.vertex))
					{
						string edgeString = RemoveEdge(L, L1);
						newGeneratedEdges.Remove(edgeString);
						VertexInfo tmp = L2;
						L2 = FindPrevVertex(L, L1, false);
						L1 = tmp;
					}
				}
				else
				{
					B = true;
				}

				if(A)
				{
					L = L1;
				}
				else if(B)
				{
					R = R1;
				}
				else if(TriangleInfo.IsAPointInsideTrianglesCircumcircle(L1.vertex, L.vertex, R.vertex, R1.vertex))
				{
					R = R1;
				}
				else
				{
					L = L1;
				}

				lowBoundEdge = new KeyValuePair<VertexInfo, VertexInfo>(L, R);
			}
		}
		return newGeneratedEdges;
	}

	KeyValuePair<VertexInfo, VertexInfo> FindHullEdge(List<VertexInfo> leftVertices, Dictionary<string, EdgeInfo> leftEdges, List<VertexInfo> rightVertices, Dictionary<string, EdgeInfo> rightEdges, bool lowerBound)
	{
		VertexInfo X;
		VertexInfo Y;
		VertexInfo Z;
		VertexInfo Z_;
		VertexInfo Z__;
		if (lowerBound) {
			X = leftVertices [leftVertices.Count - 1];
			Y = rightVertices [0];
			Z = FirstNextPointOnContour (X, leftVertices);
			Z_ = FirstNextPointOnContour (Y, rightVertices);
			Z__ = FindPrevVertex (X, Z_, false);
		} else {
			X = rightVertices [0];
			Y = leftVertices [leftVertices.Count - 1];
			Z = FirstNextPointOnContour (X, rightVertices);
			Z_ = FirstNextPointOnContour (Y, leftVertices);
			Z__ = FindPrevVertex (X, Z_, false);
		}
		while (true) {
			if(IsRightOfVector(X, Y, Z))
			{
				VertexInfo tmp = Z;
				Z = FindPrevVertex(Z, Y, false);
				Y = tmp;
			}
			else
			{
				if(IsRightOfVector(X, Y, Z__))
				{
					VertexInfo tmp = Z__;
					Z__ = FindPrevVertex(Z__, X, true);
					X = tmp;
				}
				else
				{
					return new KeyValuePair<VertexInfo, VertexInfo> (X, Y);
				}
			}
		}
	}

	VertexInfo FirstNextPointOnContour(VertexInfo vertex, List<VertexInfo> vertices)
	{
		for (int i = 0; i < vertices.Count; i++) {
			if(vertex.index == vertices[i].index)
				continue;
			if(IsNextContourPoint(vertex, vertices[i], vertices, cutPlane.normal))
				return vertices[i];
		}
		return vertex;

	}

	bool IsRightOfVector(VertexInfo vertexCenter, VertexInfo vertexEnd, VertexInfo vertex)
	{
		return Vector3.Dot (Vector3.Cross (vertexEnd.vertex - vertexCenter.vertex, vertex.vertex - vertexCenter.vertex), cutPlane.normal) > 0;
	}

	VertexInfo FindPrevVertex(VertexInfo vertexCenter, VertexInfo vertexEnd, bool isCCW)
	{
		List<string> adjancentEdges = vertexCenter.belongToEdgeIndex;

		if (!adjancentEdges.Contains (EdgeInfo.GetEdgeString (vertexCenter.index, vertexEnd.index)))
			throw new Exception("No such edge exist");
		if (adjancentEdges.Count == 1)
			return vertexEnd;
		EdgeInfo currentEdge = edges [EdgeInfo.GetEdgeString (vertexCenter.index, vertexEnd.index)];
		List<string> rightEdges = adjancentEdges.Where (e => Vector3.Dot (Vector3.Cross (vertexEnd.vertex - vertexCenter.vertex, vertices [currentEdge.GetOtherPoint (vertexCenter.index)].vertex - vertexCenter.vertex), cutPlane.normal) > 0).ToList();
		List<string> leftEdges = adjancentEdges.Except (rightEdges).ToList ();
		string nextEdgeString;
		Func<string, float> calculateDot = eString => {EdgeInfo e = edges[eString]; return Vector3.Dot((vertexEnd.vertex - vertexCenter.vertex).normalized, (vertices [e.GetOtherPoint (vertexCenter.index)].vertex - vertexCenter.vertex).normalized);};
		if (!isCCW) {
			if(rightEdges.Count>0)
			{
				nextEdgeString = rightEdges.OrderByDescending(calculateDot).First();
			}
			else
			{
				nextEdgeString = leftEdges.OrderBy(calculateDot).First();
			}
		} else {
			if(leftEdges.Count>0)
			{
				nextEdgeString = leftEdges.OrderByDescending(calculateDot).First();
			}
			else
			{
				nextEdgeString = rightEdges.OrderBy(calculateDot).First();
			}
		}

		return vertices[edges[nextEdgeString].GetOtherPoint(vertexCenter.index)];
	}
	#endregion
}

