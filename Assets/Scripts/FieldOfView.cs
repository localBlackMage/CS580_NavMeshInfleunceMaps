using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class FieldOfView : MonoBehaviour {

	Transform head;
	int layerToCastOn;
	List<PoissonDot> dots;

	// Use this for initialization
	void Start () {
		head = transform.parent;
		layerToCastOn = LayerMask.GetMask("BlockLoS");
		dots = new List<PoissonDot>();
	}
	
	// Update is called once per frame
	void Update () {
		foreach (PoissonDot pDot in dots)
		{
			Vector3 dir = pDot.transform.position - head.position;
			Ray ray = new Ray(head.position, dir.normalized);
			RaycastHit hit;
			Physics.Raycast(ray, out hit, dir.magnitude, layerToCastOn);
			if (hit.collider)
			{
				PoissonDot otherDot = hit.collider.GetComponent<PoissonDot>();
				if (otherDot == pDot)
				{
					//Debug.DrawLine(transform.parent.position, pDot.transform.position, Color.green);
					pDot.InFieldOfView();
				}
			}
		}
	}

	void OnTriggerEnter(Collider other)
	{
		PoissonDot pDot = other.GetComponent<PoissonDot>();
		if (pDot)
		{
			dots.Add(pDot);
		}
	}

	//void OnTriggerStay(Collider other)
	//{
		//PoissonDot pDot = other.GetComponent<PoissonDot>();
		//if (pDot) {
		//	Vector3 dir = other.transform.position - head.position;
		//	Ray ray = new Ray(head.position, dir.normalized);
		//	RaycastHit hit;
		//	Physics.Raycast(ray, out hit, dir.magnitude, layerToCastOn);
		//	if (hit.collider && hit.collider.Equals(other))
		//		pDot.InFieldOfView();
				
			
		//	//Debug.DrawLine(transform.parent.position, other.transform.position, Color.green);
		//}
	//}

	void OnTriggerExit(Collider other)
	{
		PoissonDot pDot = other.GetComponent<PoissonDot>();
		if (pDot)
		{
			dots.Remove(pDot);
		}
	}
}
