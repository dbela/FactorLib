
// FactorMeDlg.cpp : implementation file
//

#include "stdafx.h"
#include "FactorMe.h"
#include "FactorMeDlg.h"
#include "afxdialogex.h"
#include "factorlib.h"
#include <fstream>
#include <future>


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
	m_FactorizationTxt   = "A faktorizációs algoritmus megadja egy adott szám prímtényezős felbontását. A futási idő 30 számjegy és afelett nagy mértékben nőhet. ";
	m_DivisorFinderTxt   = "Az osztókereső algoritmus megadja egy adott szám két osztóját a kiválasztott algoritmus felhasználásával.";
	m_FermatTxt			 = "Fermat algoritmusa a megadott szám gyöke körül kezd el keresni, és akkor talál gyorsan megoldást ha a gyök körül létezik két osztója a számnak";
	m_PMinusTxt			 = "Pollard p-1 algoritmusa arra épít, hogy a megadott számnak van olyan prímosztója, hogy az attól eggyel kisebb számnak az osztói mind 10000-nél kisebbek. Ha ez teljesül nagyon gyorsan ad megoldást.";
	m_PollardRhoTxt		 = "A Pollard-Rho algoritmus a legmegbízhatóbb algoritmus 5-20 számjegyig. 20-nál nagyobb számjegyekre a futási idő miatt nem ajánlott.";
	m_QuadraticSieveTxt  = "A kvadratikus szita 20 számjegy felett is képes megoldást találni, de 20 vagy annál kisebb számjegyű számokra, viszont nem olyan hatékony. Megadhatja a paramétereket is, a B az alapul vett faktorbázisért felel, az M a szita intervallumának hosszáért. Ha az egyiket megadta, meg kell adnia a másikat is. 5-nél kisebb számokra csak paraméterekkel együtt fut le.";
	m_PrimeTestTxt       = "A prímtesztelő algoritmus megmondja a megadott számról, hogy prím-e. A maximum megadható érték: 18,446,744,073,709,551,615.";
 	m_PrimeEnumeratorTxt = "A prímsoroló algoritmus kiírja a kiválasztott fájlba a megadott számig található prímek listáját. A maximum megadható érték: 18,446,744,073,709,551,615.";
}

void CFactorMeDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Control( pDX, IDC_TAB1          , m_TabCtrl          );
	DDX_Control( pDX, IDC_RADIO2        , m_RadioButton2     );
	DDX_Control( pDX, IDC_RADIO3        , m_RadioButton3     );
	DDX_Control( pDX, IDC_RADIO4        , m_RadioButton4     );
	DDX_Control( pDX, IDC_RADIO5        , m_RadioButton5     );
	DDX_Control( pDX, IDC_STATIC1       , m_TabDetails       );
	DDX_Control( pDX, IDC_STATIC2       , m_AlgorithmDetails );
	DDX_Control( pDX, IDC_STATIC4       , m_Result           );
	DDX_Control( pDX, IDC_STATIC5       , m_ChooseFile       );
	DDX_Control( pDX, IDC_STATIC6       , m_ChooseB          );
	DDX_Control( pDX, IDC_STATIC7       , m_ChooseM          );
	DDX_Control( pDX, IDC_EDIT1         , m_EditCtrlNum      );
	DDX_Control( pDX, IDC_EDIT2         , m_EditCtrlB        );
	DDX_Control( pDX, IDC_EDIT3         , m_EditCtrlM        );
	DDX_Control( pDX, IDC_BUTTON1       , m_StartButton		 );
	DDX_Control( pDX, IDC_MFCEDITBROWSE1, m_FileBrowseCtrl   );
	DDX_Radio  ( pDX, IDC_RADIO2        , m_AlgorithmType    );
}

BEGIN_MESSAGE_MAP(CFactorMeDlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_NOTIFY(TCN_SELCHANGE, IDC_TAB1, &CFactorMeDlg::OnTcnSelchangeTab1)
	ON_BN_CLICKED(IDC_RADIO2, &CFactorMeDlg::OnBnClickedRadio2)
	ON_BN_CLICKED(IDC_RADIO3, &CFactorMeDlg::OnBnClickedRadio3)
	ON_BN_CLICKED(IDC_RADIO4, &CFactorMeDlg::OnBnClickedRadio4)
	ON_BN_CLICKED(IDC_RADIO5, &CFactorMeDlg::OnBnClickedRadio5)
	ON_BN_CLICKED(IDC_BUTTON1, &CFactorMeDlg::OnBnClickedStartButton1)
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

	m_RadioButton2    .ShowWindow( SW_HIDE );
	m_RadioButton3    .ShowWindow( SW_HIDE );
	m_RadioButton4    .ShowWindow( SW_HIDE );
	m_RadioButton5    .ShowWindow( SW_HIDE );
	m_AlgorithmDetails.ShowWindow( SW_HIDE );
	m_FileBrowseCtrl  .ShowWindow( SW_HIDE );
	m_ChooseFile      .ShowWindow( SW_HIDE );
	m_ChooseB         .ShowWindow( SW_HIDE );
	m_ChooseM         .ShowWindow( SW_HIDE );
	m_EditCtrlB       .ShowWindow( SW_HIDE );
	m_EditCtrlM       .ShowWindow( SW_HIDE );

	m_EditCtrlNum.LimitText( 60 );

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
	m_EditCtrlNum.SetWindowTextW( CString("") );
	m_Result  .SetWindowTextW( CString("") );

	int iSel = m_TabCtrl.GetCurSel();

	if( iSel == 0 )
	{
		m_EditCtrlNum.LimitText( 60 );
		m_TabDetails.SetWindowTextW( m_FactorizationTxt );
	}
	else if( iSel == 1 )
	{
		m_EditCtrlNum.LimitText( 60 );
		m_TabDetails      .SetWindowTextW( m_DivisorFinderTxt );
		m_AlgorithmDetails.SetWindowTextW( m_FermatTxt        );
	}
	else if( iSel == 2 )
	{
		m_EditCtrlNum.LimitText( 20 );
		m_TabDetails.SetWindowTextW( m_PrimeTestTxt );
	}
	else if( iSel == 3 )
	{
		m_EditCtrlNum.LimitText( 20 );
		m_TabDetails.SetWindowTextW( m_PrimeEnumeratorTxt );
	}

	if( iSel == 1 )
	{
		m_RadioButton2    .ShowWindow( SW_SHOW );
		m_RadioButton3    .ShowWindow( SW_SHOW );
		m_RadioButton4    .ShowWindow( SW_SHOW );
		m_RadioButton5    .ShowWindow( SW_SHOW );
		m_AlgorithmDetails.ShowWindow( SW_SHOW );
		UpdateData( true );
		if( m_AlgorithmType == 3 )
		{
			m_ChooseB         .ShowWindow( SW_SHOW );
			m_ChooseM         .ShowWindow( SW_SHOW );
			m_EditCtrlB       .ShowWindow( SW_SHOW );
			m_EditCtrlM       .ShowWindow( SW_SHOW );
		}
		else
		{
			m_ChooseB         .ShowWindow( SW_HIDE );
			m_ChooseM         .ShowWindow( SW_HIDE );
			m_EditCtrlB       .ShowWindow( SW_HIDE );
			m_EditCtrlM       .ShowWindow( SW_HIDE );
		}
	}
	else
	{
		m_RadioButton2    .ShowWindow( SW_HIDE );
		m_RadioButton3    .ShowWindow( SW_HIDE );
		m_RadioButton4    .ShowWindow( SW_HIDE );
		m_RadioButton5    .ShowWindow( SW_HIDE );
		m_AlgorithmDetails.ShowWindow( SW_HIDE );
		m_ChooseB         .ShowWindow( SW_HIDE );
		m_ChooseM         .ShowWindow( SW_HIDE );
		m_EditCtrlB       .ShowWindow( SW_HIDE );
		m_EditCtrlM       .ShowWindow( SW_HIDE );
	}

	if( iSel == 3 )
	{
		m_FileBrowseCtrl.ShowWindow( SW_SHOW );
		m_FileBrowseCtrl.EnableFileBrowseButton(TEXT("txt"), TEXT("Text Files|*.txt||"));
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

void CFactorMeDlg::OnBnClickedRadio2()
{
	m_AlgorithmDetails.SetWindowTextW( m_FermatTxt );
	m_Result          .SetWindowTextW( CString("") );
	m_ChooseB         .ShowWindow( SW_HIDE );
	m_ChooseM         .ShowWindow( SW_HIDE );
	m_EditCtrlB       .ShowWindow( SW_HIDE );
	m_EditCtrlM       .ShowWindow( SW_HIDE );
}

void CFactorMeDlg::OnBnClickedRadio3()
{
	m_AlgorithmDetails.SetWindowTextW( m_PMinusTxt );
	m_Result          .SetWindowTextW( CString("") );
	m_ChooseB         .ShowWindow( SW_HIDE );
	m_ChooseM         .ShowWindow( SW_HIDE );
	m_EditCtrlB       .ShowWindow( SW_HIDE );
	m_EditCtrlM       .ShowWindow( SW_HIDE );
}

void CFactorMeDlg::OnBnClickedRadio4()
{
	m_AlgorithmDetails.SetWindowTextW( m_PollardRhoTxt );
	m_Result          .SetWindowTextW( CString("") );
	m_ChooseB         .ShowWindow( SW_HIDE );
	m_ChooseM         .ShowWindow( SW_HIDE );
	m_EditCtrlB       .ShowWindow( SW_HIDE );
	m_EditCtrlM       .ShowWindow( SW_HIDE );
}

void CFactorMeDlg::OnBnClickedRadio5()
{
	m_AlgorithmDetails.SetWindowTextW( m_QuadraticSieveTxt );
	m_Result          .SetWindowTextW( CString("") );
	m_ChooseB         .ShowWindow( SW_SHOW );
	m_ChooseM         .ShowWindow( SW_SHOW );
	m_EditCtrlB       .ShowWindow( SW_SHOW );
	m_EditCtrlM       .ShowWindow( SW_SHOW );
}


void CFactorMeDlg::OnBnClickedStartButton1()
{
	THREADSTRUCT *_param = new THREADSTRUCT;
    _param->_this = this;
	AfxBeginThread( DoAlgorithms, (LPVOID)_param );

}

UINT __cdecl CFactorMeDlg::DoAlgorithms( LPVOID pParam )
{
	THREADSTRUCT*    ts = (THREADSTRUCT*)pParam;
	CFactorMeDlg* Dlg = ts->_this;
	
	Dlg->EnableControls( false );
	Dlg->BeginWaitCursor();

	int iSel = Dlg->m_TabCtrl.GetCurSel();

	CString cStrNum;
	Dlg->m_EditCtrlNum.GetWindowText( cStrNum );

	CT2CA conNum( cStrNum );

	std::string strNum( conNum );

	std::string ResultMessage;

	if( strNum == "" )
	{
		ResultMessage = "Kérem adjon meg egy számot! ";
	}
	else if( strNum == "0" || strNum == "1" )
	{
		ResultMessage = strNum + " = " + strNum;
	}
	else if( strNum.find_first_not_of( "0123456789" ) != std::string::npos )
	{
		ResultMessage = "Hibás bemenet. A megadott érték negatív, vagy nem egész szám. ";
	}
	else if( iSel == 0 )
	{
		ResultMessage = FactorLib::FactorLib::RunFactorize( strNum );
	}
	else if( iSel == 1 )
	{
		Dlg->UpdateData( true );
		switch( Dlg->m_AlgorithmType )
		{
			case 0:
				ResultMessage = FactorLib::FactorLib::RunFermat( strNum );
				break;
			case 1:
				ResultMessage = FactorLib::FactorLib::RunPMinus( strNum );
				break;
			case 2:
				ResultMessage = FactorLib::FactorLib::RunPollardRho( strNum );
				break;
			case 3:
				

				CString cStrB;
				Dlg->m_EditCtrlB.GetWindowText( cStrB );
				CT2CA conB( cStrB );
				std::string strB ( conB );

				if( strB.find_first_not_of( "0123456789" ) != std::string::npos )
				{
					ResultMessage = "Hibás bemenet. A B paraméter negatív, vagy nem egész szám. ";
					break;
				}
				
				CString cStrM;
				Dlg->m_EditCtrlM.GetWindowText( cStrM );
				CT2CA conM( cStrM );
				std::string strM ( conM );

				
				if( strM.find_first_not_of( "0123456789" ) != std::string::npos )
				{
					ResultMessage = "Hibás bemenet. Az M paraméter negatív, vagy nem egész szám. ";
					break;
				}

				if( ( strB == "" && strM != "" ) || ( strB != "" && strM == "" ) )
				{
					ResultMessage = "Kérem adja meg mindkét paramétert! ";
					break;
				}

				if( strNum.length() < 5 && strB == "" )
				{
					ResultMessage = "A szám túl kicsi";
					break;
				}

				

				ResultMessage = FactorLib::FactorLib::RunQuadraticSieve( strNum, strB, strM );
				break;
		}
	}
	else if( iSel == 2 )
	{
		ResultMessage = FactorLib::FactorLib::RunPrimeTest( strNum );
	}
	else if( iSel == 3 )
	{
		CString cStrFileName;
		Dlg->m_FileBrowseCtrl.GetWindowText( cStrFileName );

		if( cStrFileName == "" )
		{
			ResultMessage = "Kérem válasszon ki egy fájlt!";
		}
		else
		{
			CT2CA conFileName( cStrFileName );

			std::string strFileName( conFileName ); 

			std::future<std::vector<unsigned long long>> result = std::async(std::launch::async, FactorLib::FactorLib::RunSieveOfE, strNum );
			std::vector<unsigned long long> primes = result.get();
		
			if( primes.size() > 0)
			{
				std::ofstream file;
				file.open( strFileName );
				if( file.is_open() )
				{
					for( int i = 0; i < primes.size(); ++i )
					{
						file << primes[i] << " ";
						if( i % 5 == 0 && i != 0 )
						{
							file << "\n";
						}
					}
					ResultMessage = "A fájl kiírás sikeres.";
				}
				else
				{
					ResultMessage = "A megadott fájl nem elérhető.";
				}
			}
			else
			{
				ResultMessage = "A szám nem a megadott értékeken belül van.";
			}
		}
	}

	CString resMsg( ResultMessage.c_str() );
	Dlg->m_Result.SetWindowTextW( resMsg );

	Dlg->EnableControls( true );
	Dlg->EndWaitCursor();

	return 1;
}

void CFactorMeDlg::EnableControls( bool bTo )
{
	m_TabCtrl         .EnableWindow(bTo);
	m_RadioButton2    .EnableWindow(bTo);
	m_RadioButton3    .EnableWindow(bTo);
	m_RadioButton4    .EnableWindow(bTo);
	m_RadioButton5    .EnableWindow(bTo);
	m_FileBrowseCtrl  .EnableWindow(bTo);
	m_StartButton     .EnableWindow(bTo);
	m_EditCtrlNum     .EnableWindow(bTo);

}

