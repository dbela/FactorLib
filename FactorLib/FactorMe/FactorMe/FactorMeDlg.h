
// FactorMeDlg.h : header file
//

#pragma once



// CFactorMeDlg dialog
class CFactorMeDlg : public CDialogEx
{
// Construction
public:
	CFactorMeDlg(CWnd* pParent = nullptr);	// standard constructor

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_FACTORME_DIALOG };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support


// Implementation
protected:
	HICON m_hIcon;
	CTabCtrl           m_TabCtrl;
	int                m_AlgorithmType;
	CButton            m_RadioButton2;
	CButton            m_RadioButton3;
	CButton            m_RadioButton4;
	CButton            m_RadioButton5;
	CStatic            m_TabDetails;
	CStatic            m_AlgorithmDetails;
	CStatic            m_Result;
	CStatic            m_ChooseFile;
	CStatic            m_ChooseB;
	CStatic            m_ChooseM;
	CMFCEditBrowseCtrl m_FileBrowseCtrl;
	CButton            m_StartButton;
	CEdit              m_EditCtrlNum;
	CEdit              m_EditCtrlB;
	CEdit              m_EditCtrlM;

	CString            m_FactorizationTxt;
	CString            m_DivisorFinderTxt;
	CString            m_PrimeTestTxt;
	CString            m_PrimeEnumeratorTxt;
	CString            m_FermatTxt;
	CString            m_PMinusTxt;
	CString            m_PollardRhoTxt;
	CString            m_QuadraticSieveTxt;

	static UINT DoAlgorithms( LPVOID pParam );
	void EnableControls( bool bTo );

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnTcnSelchangeTab1(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnBnClickedRadio4();
	afx_msg void OnBnClickedRadio3();
	afx_msg void OnBnClickedRadio2();
	afx_msg void OnBnClickedRadio5();
	afx_msg void OnBnClickedStartButton1();
};

typedef struct THREADSTRUCT
{
    CFactorMeDlg*    _this;
        //you can add here other parameters you might be interested on
} THREADSTRUCT;
