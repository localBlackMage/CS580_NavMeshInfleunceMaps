using UnityEngine;
using System.Collections;

public class NavCam : MonoBehaviour
{
	public float cameraSensitivity = 90;
	public float climbSpeed = 4;
	public float normalMoveSpeed = 10;

	private float rotationX = 0.0f;
	private float rotationY = 0.0f;

	void Start()
	{
		Cursor.lockState = CursorLockMode.Confined;
	}

	void Update()
	{
		rotationX += Input.GetAxis("Mouse X") * cameraSensitivity * Time.deltaTime;
		rotationY += Input.GetAxis("Mouse Y") * cameraSensitivity * Time.deltaTime;
		rotationY = Mathf.Clamp(rotationY, -90, 90);

		transform.localRotation = Quaternion.AngleAxis(rotationX, Vector3.up);
		transform.localRotation *= Quaternion.AngleAxis(rotationY, Vector3.left);

		transform.position += transform.forward * normalMoveSpeed * Input.GetAxis("Vertical") * Time.deltaTime;
		transform.position += transform.right * normalMoveSpeed * Input.GetAxis("Horizontal") * Time.deltaTime;

		if (Input.GetKey(KeyCode.E)) { transform.position += transform.up * climbSpeed * Time.deltaTime; }
		if (Input.GetKey(KeyCode.Q)) { transform.position -= transform.up * climbSpeed * Time.deltaTime; }

		if (Input.GetKeyDown(KeyCode.End))
		{
			Cursor.lockState = Cursor.lockState == CursorLockMode.Confined ? CursorLockMode.Locked : CursorLockMode.Confined;
		}
	}
}
