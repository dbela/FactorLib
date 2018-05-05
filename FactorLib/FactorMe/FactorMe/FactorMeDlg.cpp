
// FactorMeDlg.cpp : implementation file
//

#include "stdafx.h"
#include "FactorMe.h"
#include "FactorMeDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CAboutDlg dialog used for App About

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_ABOUTBOX };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(IDD_ABOUTBOX)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// CFactorMeDlg dialog



CFactorMeDlg::CFactorMeDlg(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_FACTORME_DIALOG, pParent), m_AlgorithmType(0)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
	m_FactorizationTxt   = "A faktorizációs algoritmus megadja egy adott szám prímtényezős felbontását. Az algoritmus csak 40 számjegyű vagy azalatti számokra működik.";
	m_DivisorFinderTxt   = "Az osztókereső algoritmus megadja egy adott szám két osztóját a kiválasztott algoritmus felhasználásával.";
	m_TrialDivTxt        = "A próbaosztás a legegyszerűbb osztókereső algoritmus, viszont minél nagyobb a szám, annál tovább tarthat az osztó keresése. 10000-nél nagyobb számokra már nem ajánlatos";
	m_FermatTxt			 = "Fermat algoritmusa a megadott szám gyöke körül kell keresni az osztókat, így akkor talál gyorsan megoldást ha a gyök körül létezik osztó";
	m_PMinusTxt			 = "Pollard p-1 algoritmusa arra épít, hogy a megadott számtól eggyel kisebb számnak az osztói mind 10000-nél kisebbek. Ha ez nem teljesül, nem kapunk megoldást.";
	m_PollardRhoTxt		 = "A Pollard-Rho algoritmus a legmegbízhatóbb algoritmus 5-20 számjegyig. 20-nál nagyobb számjegyekre a futási idő miatt nem ajánlott.";
	m_QuadraticSieveTxt  = "A kvadratikus szita az igazi nagyágyú az itt látható algoritmusok közül, hiszen 20-40 számjegyig is képes aránylag gyorsan megoldást találni.";
	m_PrimeTestTxt       = "A prímtesztelő algoritmus megmondja a megadott számról, hogy prím-e. A maximum megadható érték: 18,446,744,073,709,551,615.";
 	m_PrimeEnumeratorTxt = "A prímsoroló algoritmus kiírja a megadott fájlba a megadott számig található prímek listáját. A maximum megadható érték: 18,446,744,073,709,551,615.";
}

void CFactorMeDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Control( pDX, IDC_TAB1          , m_TabCtrl          );
	DDX_Control( pDX, IDC_RADIO1        , m_RadioButton1     );
	DDX_Control( pDX, IDC_RADIO2        , m_RadioButton2     );
	DDX_Control( pDX, IDC_RADIO3        , m_RadioButton3     );
	DDX_Control( pDX, IDC_RADIO4        , m_RadioButton4     );
	DDX_Control( pDX, IDC_RADIO5        , m_RadioButton5     );
	DDX_Control( pDX, IDC_STATIC1       , m_TabDetails       );
	DDX_Control( pDX, IDC_STATIC2       , m_AlgorithmDetails );
	DDX_Control( pDX, IDC_STATIC4       , m_Result           );
	DDX_Control( pDX, IDC_STATIC5       , m_ChooseFile       );
	DDX_Control( pDX, IDC_EDIT1         , m_EditCtrl         );
	DDX_Control( pDX, IDC_BUTTON1       , m_StartButton		 );
	DDX_Control( pDX, IDC_MFCEDITBROWSE1, m_FileBrowseCtrl   );
	DDX_Radio  ( pDX, IDC_RADIO1        , m_AlgorithmType    );
}

BEGIN_MESSAGE_MAP(CFactorMeDlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_NOTIFY(TCN_SELCHANGE, IDC_TAB1, &CFactorMeDlg::OnTcnSelchangeTab1)
	ON_BN_CLICKED(IDC_RADIO1, &CFactorMeDlg::OnBnClickedRadio1)
	ON_BN_CLICKED(IDC_RADIO2, &CFactorMeDlg::OnBnClickedRadio2)
	ON_BN_CLICKED(IDC_RADIO3, &CFactorMeDlg::OnBnClickedRadio3)
	ON_BN_CLICKED(IDC_RADIO4, &CFactorMeDlg::OnBnClickedRadio4)
	ON_BN_CLICKED(IDC_RADIO5, &CFactorMeDlg::OnBnClickedRadio5)
END_MESSAGE_MAP()


// CFactorMeDlg message handlers

BOOL CFactorMeDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != nullptr)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	m_TabCtrl.InsertItem( 0, _T( "Faktorizáció" ) );
	m_TabCtrl.InsertItem( 1, _T( "Osztókereső" ) );
	m_TabCtrl.InsertItem( 2, _T( "Prímtesztelő" ) );
	m_TabCtrl.InsertItem( 3, _T( "Prímsoroló" ) );

	m_RadioButton1    .ShowWindow( SW_HIDE );
	m_RadioButton2    .ShowWindow( SW_HIDE );
	m_RadioButton3    .ShowWindow( SW_HIDE );
	m_RadioButton4    .ShowWindow( SW_HIDE );
	m_RadioButton5    .ShowWindow( SW_HIDE );
	m_AlgorithmDetails.ShowWindow( SW_HIDE );
	m_FileBrowseCtrl  .ShowWindow( SW_HIDE );
	m_ChooseFile      .ShowWindow( SW_HIDE );

	m_TabDetails.SetWindowTextW( m_FactorizationTxt );
	m_Result    .SetWindowTextW( CString("") );

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	// TODO: Add extra initialization here

	

	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CFactorMeDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CFactorMeDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CFactorMeDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}



void CFactorMeDlg::OnTcnSelchangeTab1(NMHDR *pNMHDR, LRESULT *pResult)
{
	m_Result.SetWindowTextW( CString("") );

	int iSel = m_TabCtrl.GetCurSel();

	if( iSel == 0 )
	{
		m_TabDetails.SetWindowTextW( m_FactorizationTxt );
	}
	else if( iSel == 1 )
	{
		m_TabDetails      .SetWindowTextW( m_DivisorFinderTxt );
		m_AlgorithmDetails.SetWindowTextW( m_TrialDivTxt      );
	}
	else if( iSel == 2 )
	{
		m_TabDetails.SetWindowTextW( m_PrimeTestTxt );
	}
	else
	{
		m_TabDetails.SetWindowTextW( m_PrimeEnumeratorTxt );
	}

	if( iSel == 1 )
	{
		m_RadioButton1    .ShowWindow( SW_SHOW );
		m_RadioButton2    .ShowWindow( SW_SHOW );
		m_RadioButton3    .ShowWindow( SW_SHOW );
		m_RadioButton4    .ShowWindow( SW_SHOW );
		m_RadioButton5    .ShowWindow( SW_SHOW );
		m_AlgorithmDetails.ShowWindow( SW_SHOW );
	}
	else
	{
		m_RadioButton1    .ShowWindow( SW_HIDE );
		m_RadioButton2    .ShowWindow( SW_HIDE );
		m_RadioButton3    .ShowWindow( SW_HIDE );
		m_RadioButton4    .ShowWindow( SW_HIDE );
		m_RadioButton5    .ShowWindow( SW_HIDE );
		m_AlgorithmDetails.ShowWindow( SW_HIDE );
	}

	if( iSel == 3 )
	{
		m_FileBrowseCtrl.ShowWindow( SW_SHOW );
		m_ChooseFile    .ShowWindow( SW_SHOW );
	}
	else 
	{
		m_FileBrowseCtrl.ShowWindow( SW_HIDE );
		m_ChooseFile    .ShowWindow( SW_HIDE );
	}
	// TODO: Add your control notification handler code here
	*pResult = 0;
}

void CFactorMeDlg::OnBnClickedRadio1()
{
	m_AlgorithmDetails.SetWindowTextW( m_TrialDivTxt );
}

void CFactorMeDlg::OnBnClickedRadio2()
{
	m_AlgorithmDetails.SetWindowTextW( m_FermatTxt );
}

void CFactorMeDlg::OnBnClickedRadio3()
{
	m_AlgorithmDetails.SetWindowTextW( m_PMinusTxt );
}

void CFactorMeDlg::OnBnClickedRadio4()
{
	m_AlgorithmDetails.SetWindowTextW( m_PollardRhoTxt );
}

void CFactorMeDlg::OnBnClickedRadio5()
{
	m_AlgorithmDetails.SetWindowTextW( m_QuadraticSieveTxt );
}
